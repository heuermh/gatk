package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Reference implementation of multi-allelic EXACT model.  Extremely slow for many alternate alleles.
 */
public final class ReferenceDiploidExactAFCalculator extends DiploidExactAFCalculator {

    private static final double LOG_OF_2 = MathUtils.LogCache.get(2);

    protected ReferenceDiploidExactAFCalculator() {
    }

    @Override
    protected AFCalculationResult computeLogPNonRef(final VariantContext vc, final int defaultPloidy,
                                                    final double[] logAlleleFrequencyPriors, final StateTracker stateTracker) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(logAlleleFrequencyPriors, "logAlleleFrequencyPriors is null");
        Utils.nonNull(stateTracker, "stateTracker is null");
        final int numAlternateAlleles = vc.getNAlleles() - 1;

        final List<double[]> genotypeLikelihoods = getGLs(vc.getGenotypes(), true);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        // queue of AC conformations to process
        final Deque<ExactACset> ACqueue = new LinkedList<>();

        // mapping of ExactACset indexes to the objects
        final Map<ExactACcounts, ExactACset> indexesToACset = new HashMap<>(numChr+1);

        // add AC=0 to the queue
        final int[] zeroCounts = new int[numAlternateAlleles];
        final ExactACset zeroSet = new ExactACset(numSamples+1, new ExactACcounts(zeroCounts));
        ACqueue.add(zeroSet);
        indexesToACset.put(zeroSet.getACcounts(), zeroSet);

        while ( !ACqueue.isEmpty() ) {

            // compute logLikelihoods
            final ExactACset set = ACqueue.remove();

            calculateAlleleCountConformation(set, genotypeLikelihoods, numChr, ACqueue, indexesToACset, logAlleleFrequencyPriors,stateTracker);

            // clean up memory
            indexesToACset.remove(set.getACcounts());
        }

        return getResultFromFinalState(vc, logAlleleFrequencyPriors, stateTracker);
    }


    private double calculateAlleleCountConformation(final ExactACset set,
                                                    final List<double[]> genotypeLikelihoods,
                                                    final int numChr,
                                                    final Deque<ExactACset> ACqueue,
                                                    final Map<ExactACcounts, ExactACset> indexesToACset,
                                                    final double[] logAlleleFrequencyPriors,
                                                    final StateTracker stateTracker) {

        // compute the logLikelihoods
        computeLofK(set, genotypeLikelihoods, logAlleleFrequencyPriors, stateTracker);

        final double logLofK = set.getLogLikelihoods(set.size()-1);

        // can we abort early because the logLikelihoods are so small?
        if ( stateTracker.abort(logLofK, set.getACcounts(), true, false) ) {
            return logLofK;
        }

        // iterate over higher frequencies if possible
        final int ACwiggle = numChr - set.getACsum();
        if ( ACwiggle == 0 ){ // all alternate alleles already sum to 2N so we cannot possibly go to higher frequencies
            return logLofK;
        }

        final int numAltAlleles = set.getACcounts().getCounts().length;

        // add conformations for the k+1 case
        for ( int allele = 0; allele < numAltAlleles; allele++ ) {
            final int[] ACcountsClone = set.getACcounts().getCounts().clone();
            ACcountsClone[allele]++;
            // to get to this conformation, a sample would need to be AB (remember that ref=0)
            final int PLindex = GenotypeLikelihoods.calculatePLindex(0, allele + 1);
            updateACset(ACcountsClone, numChr, set, PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
        }

        // add conformations for the k+2 case if it makes sense; note that the 2 new alleles may be the same or different
        if ( ACwiggle > 1 ) {
            final List<DependentSet> differentAlleles = new ArrayList<>(numAltAlleles * numAltAlleles);
            final List<DependentSet> sameAlleles = new ArrayList<>(numAltAlleles);

            for ( int allele_i = 0; allele_i < numAltAlleles; allele_i++ ) {
                for ( int allele_j = allele_i; allele_j < numAltAlleles; allele_j++ ) {
                    final int[] ACcountsClone = set.getACcounts().getCounts().clone();
                    ACcountsClone[allele_i]++;
                    ACcountsClone[allele_j]++;

                    // to get to this conformation, a sample would need to be BB or BC (remember that ref=0, so add one to the index)
                    final int PLindex = GenotypeLikelihoods.calculatePLindex(allele_i + 1, allele_j + 1);
                    if ( allele_i == allele_j ) {
                        sameAlleles.add(new DependentSet(ACcountsClone, PLindex));
                    } else {
                        differentAlleles.add(new DependentSet(ACcountsClone, PLindex));
                    }
                }
            }

            // IMPORTANT: we must first add the cases where the 2 new alleles are different so that the queue maintains its ordering
            for ( final DependentSet dependent : differentAlleles ) {
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
            }
            for ( final DependentSet dependent : sameAlleles ) {
                updateACset(dependent.ACcounts, numChr, set, dependent.PLindex, ACqueue, indexesToACset, genotypeLikelihoods);
            }
        }

        return logLofK;
    }

    private static void computeLofK(final ExactACset set,
                                    final List<double[]> genotypeLikelihoods,
                                    final double[] logAlleleFrequencyPriors, final StateTracker stateTracker) {

        set.setLogLikelihoods(0, 0.0); // the zero case
        final int totalK = set.getACsum();

        // special case for k = 0 over all k
        if ( totalK == 0 ) {
            for ( int j = 1; j < set.size(); j++ ) {
                set.setLogLikelihoods(j, set.getLogLikelihoods(j - 1) + genotypeLikelihoods.get(j)[HOM_REF_INDEX]);
            }

            final double logLof0 = set.getLogLikelihoods(set.size()-1);
            stateTracker.setLogLikelihoodOfAFzero(logLof0);
            stateTracker.setLogPosteriorOfAFzero(logLof0 + logAlleleFrequencyPriors[0]);
            return;
        }

        // if we got here, then k > 0 for at least one k.
        // the non-AA possible conformations were already dealt with by pushes from dependent sets;
        // now deal with the AA case (which depends on previous cells in this column) and then update the L(j,k) value
        for ( int j = 1; j < set.size(); j++ ) {
            if ( totalK < 2*j-1 ) {
                final double[] gl = genotypeLikelihoods.get(j);
                final double conformationValue = MathUtils.LogCache.get(2*j-totalK) + MathUtils.LogCache.get(2*j-totalK-1) + set.getLogLikelihoods(j-1) + gl[HOM_REF_INDEX];
                set.setLogLikelihoods(j, MathUtils.approximateLogSumLog(set.getLogLikelihoods(j), conformationValue));
            }

            final double logDenominator = MathUtils.LogCache.get(2*j) + MathUtils.LogCache.get(2*j-1);
            set.setLogLikelihoods(j, set.getLogLikelihoods(j) - logDenominator);
        }

        double logLofK = set.getLogLikelihoods(set.size()-1);

        // update the MLE if necessary
        stateTracker.updateMLEifNeeded(logLofK, set.getACcounts().getCounts());

        // apply the priors over each alternate allele
        for ( final int ACcount : set.getACcounts().getCounts() ) {
            if ( ACcount > 0 ) {
                logLofK += logAlleleFrequencyPriors[ACcount];
            }
        }

        stateTracker.updateMAPifNeeded(logLofK, set.getACcounts().getCounts());
    }

    private static final class DependentSet {
        public final int[] ACcounts;
        public final int PLindex;

        DependentSet(final int[] ACcounts, final int PLindex) {
            this.ACcounts = ACcounts;
            this.PLindex = PLindex;
        }
    }


    // adds the ExactACset represented by the ACcounts to the ACqueue if not already there (creating it if needed) and
    // also pushes its value to the given callingSetIndex.
    private static void updateACset(final int[] newSetCounts,
                                    final int numChr,
                                    final ExactACset dependentSet,
                                    final int PLsetIndex,
                                    final Queue<ExactACset> ACqueue,
                                    final Map<ExactACcounts, ExactACset> indexesToACset,
                                    final List<double[]> genotypeLikelihoods) {
        final ExactACcounts index = new ExactACcounts(newSetCounts);
        if ( !indexesToACset.containsKey(index) ) {
            final ExactACset set = new ExactACset(numChr/2 +1, index);
            indexesToACset.put(index, set);
            ACqueue.add(set);
        }

        // push data from the dependency to the new set
        pushData(indexesToACset.get(index), dependentSet, PLsetIndex, genotypeLikelihoods);
    }

    private static void pushData(final ExactACset targetSet,
                                 final ExactACset dependentSet,
                                 final int PLsetIndex,
                                 final List<double[]> genotypeLikelihoods) {
        final int totalK = targetSet.getACsum();

        for ( int j = 1; j < targetSet.size(); j++ ) {
            if ( totalK <= 2*j ) { // skip impossible conformations
                final double[] gl = genotypeLikelihoods.get(j);
                final double conformationValue =
                        determineCoefficient(PLsetIndex, j, targetSet.getACcounts().getCounts(), totalK) + dependentSet.getLogLikelihoods(j-1) + gl[PLsetIndex];
                targetSet.setLogLikelihoods(j, MathUtils.approximateLogSumLog(targetSet.getLogLikelihoods(j), conformationValue));
            }
        }
    }

    private static double determineCoefficient(final int PLindex, final int j, final int[] ACcounts, final int totalK) {
        // the closed form representation generalized for multiple alleles is as follows:
        // AA: (2j - totalK) * (2j - totalK - 1)
        // AB: 2k_b * (2j - totalK)
        // AC: 2k_c * (2j - totalK)
        // BB: k_b * (k_b - 1)
        // BC: 2 * k_b * k_c
        // CC: k_c * (k_c - 1)

        // find the 2 alleles that are represented by this PL index
        final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindex);

        // *** note that throughout this method we subtract one from the alleleIndex because ACcounts ***
        // *** doesn't consider the reference allele whereas the GenotypeLikelihoods PL cache does.   ***

        // the AX het case
        if ( alleles.alleleIndex1 == 0 ) {
            return MathUtils.LogCache.get(2 * ACcounts[alleles.alleleIndex2 - 1]) + MathUtils.LogCache.get(2 * j - totalK);
        }

        final int k_i = ACcounts[alleles.alleleIndex1-1];

        // the hom var case (e.g. BB, CC, DD)
        final double coeff;
        if ( alleles.alleleIndex1 == alleles.alleleIndex2 ) {
            coeff = MathUtils.LogCache.get(k_i) + MathUtils.LogCache.get(k_i - 1);
        } else {        // the het non-ref case (e.g. BC, BD, CD)
            final int k_j = ACcounts[alleles.alleleIndex2-1];
            coeff = LOG_OF_2 + MathUtils.LogCache.get(k_i) + MathUtils.LogCache.get(k_j);
        }

        return coeff;
    }
}

package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Arrays;

public final class BaseRecalibratorSparkFn {

    public static RecalibrationReport apply( final JavaPairRDD<GATKRead, ReadContextData> readsWithContext, final SAMFileHeader header, final SAMSequenceDictionary referenceDictionary, final RecalibrationArgumentCollection recalArgs ) {
        final ReadFilter filter = BaseRecalibrator.getStandardBQSRReadFilter(header);
        final JavaPairRDD<GATKRead, ReadContextData> filtered = readsWithContext.filter(readWithContext -> filter.apply(readWithContext._1()));

        JavaRDD<RecalibrationTables> unmergedTables = filtered.mapPartitions(readWithContextIterator -> {
            final BaseRecalibrationEngine bqsr = new BaseRecalibrationEngine(recalArgs, header);
            bqsr.logCovariatesUsed();

            while ( readWithContextIterator.hasNext() ) {
                final Tuple2<GATKRead, ReadContextData> readWithData = readWithContextIterator.next();
                Iterable<Variant> variants = readWithData._2().getOverlappingVariants();
                final ReferenceBases refBases = readWithData._2().getOverlappingReferenceBases();
                ReferenceDataSource refDS = new ReferenceMemorySource(refBases, referenceDictionary);

                bqsr.processRead(readWithData._1(), refDS, variants);
            }
            // Need to wrap in ArrayList due to our current inability to serialize the return value of Arrays.asList() directly
            return new ArrayList<>(Arrays.asList(bqsr.getRecalibrationTables()));
        });

        final RecalibrationTables emptyRecalibrationTable = new RecalibrationTables(new StandardCovariateList(recalArgs, header));
        final RecalibrationTables combinedTables = unmergedTables.treeAggregate(emptyRecalibrationTable,
                RecalibrationTables::inPlaceCombine,
                RecalibrationTables::inPlaceCombine,
                Math.max(1, (int)(Math.log(unmergedTables.partitions().size()) / Math.log(2))));

        BaseRecalibrationEngine.finalizeRecalibrationTables(combinedTables);

        final QuantizationInfo quantizationInfo = new QuantizationInfo(combinedTables, recalArgs.QUANTIZING_LEVELS);

        final StandardCovariateList covariates = new StandardCovariateList(recalArgs, header);
        final RecalibrationArgumentCollection RAC = recalArgs;
        final boolean sortByCols = recalArgs.SORT_BY_ALL_COLUMNS;
        final RecalibrationTables recalibrationTables = combinedTables;
        final GATKReport report = RecalUtils.createRecalibrationGATKReport(RAC.generateReportTable(covariates.covariateNames()), quantizationInfo.generateReportTable(sortByCols), RecalUtils.generateReportTables(recalibrationTables, covariates, sortByCols));
        return new RecalibrationReport(report);
    }
}

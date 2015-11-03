package org.broadinstitute.hellbender.utils.test;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.SamComparison;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Collection of utilities for making common assertions about SAM files for unit testing purposes.
 */
public final class SamAssertionUtils {

    private static SamReader getReader(final File sam, final ValidationStringency validationStringency, final File reference) {
        return SamReaderFactory.makeDefault().validationStringency(validationStringency).referenceSequence(reference).open(sam);
    }

    /**
     *  causes an exception if the given sam files aren't equal
     *  @param reference is allowed to be null
     */
    public static void assertSamsEqual(final File sam1, final File sam2, final ValidationStringency validationStringency, final File reference) throws IOException {
        final String equalStringent = samsEqualStringent(sam1, sam2, validationStringency, reference);
        Assert.assertNull(equalStringent, "SAM file " + sam1.getPath() + " differs from expected output:" + sam2.getPath() + " " + equalStringent);
    }

    /**
     * causes an exception if the given sam files aren't equal
     */
    public static void assertSamsEqual(final File sam1, final File sam2, final ValidationStringency validationStringency) throws IOException {
        assertSamsEqual(sam1, sam2, validationStringency, null);
    }

    /**
     * causes an exception if the given sam files aren't equal
     * @param reference is allowed to be null
     * the default ValidationStringency value for this method is DEFAULT_STRINGENCY
     */
    public static void assertSamsEqual(final File sam1, final File sam2, final File reference) throws IOException {
        assertSamsEqual(sam1, sam2, ValidationStringency.DEFAULT_STRINGENCY, reference);
    }

    /**
     * causes an exception if the given sam files aren't equal
     * the default ValidationStringency value for this method is DEFAULT_STRINGENCY
     */
    public static void assertSamsEqual(final File sam1, final File sam2) throws IOException {
        assertSamsEqual(sam1, sam2, ValidationStringency.DEFAULT_STRINGENCY, null);
    }

    /**
     * causes an exception if the given sam isn't valid
     * @param reference is allowed to be null
     */
    public static void assertSamValid(final File sam, final ValidationStringency validationStringency, final File reference) throws IOException {
        try (final SamReader samReader = getReader(sam, validationStringency, reference)) {
            final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);
            validator.setIgnoreWarnings(true);
            validator.setVerbose(true, 1000);
            validator.setErrorsToIgnore(Arrays.asList(SAMValidationError.Type.MISSING_READ_GROUP));
            final boolean validated = validator.validateSamFileVerbose(samReader, null);
            Assert.assertTrue(validated, "SAM file validation failed");
        }
    }

    /**
     * causes an exception if the given sam isn't valid
     */
    public static void assertSamValid(final File sam, final ValidationStringency validationStringency) throws IOException {
        assertSamValid(sam, validationStringency, null);
    }

    /**
     * causes an exception if the given sam isn't valid
     * @param reference is allowed to be null
     * the default ValidationStringency value for this method is LENIENT
     */
    public static void assertSamValid(final File sam, final File reference) throws IOException {
        assertSamValid(sam, ValidationStringency.LENIENT, reference);
    }

    /**
     * causes an exception if the given sam isn't valid
     * the default ValidationStringency value for this method is LENIENT
     */
    public static void assertSamValid(final File sam) throws IOException {
        assertSamValid(sam, ValidationStringency.LENIENT, null);
    }

    /**
     * Compares SAM/BAM files in a stringent way but not by byte identity (allow reorder of attributes)
     * Comparing by MD5s is too strict and comparing by SamComparison is too lenient. So we need this method.
     * Returns null if equal or message string if not equal.
     */
    public static String samsEqualStringent(final File sam1, final File sam2, final ValidationStringency validation, final File reference) throws IOException {
        //Three phases:
        //1:  if MD5s are same, then we're equal. Quick check
        final String fileMD5_1 = Utils.calculateFileMD5(sam1);
        final String fileMD5_2 = Utils.calculateFileMD5(sam2);
        if (fileMD5_1.equals(fileMD5_2)){
            return null;
        }

        //2:  SamComparison says no, then we're not equal
        try(final SamReader reader1 = getReader(sam1, validation, reference);
            final SamReader reader2 = getReader(sam2, validation, reference)){

            final SamComparison comparison = new SamComparison(reader1, reader2);
            if (! comparison.areHeadersEqual()){  //For headers, we want identity (at least for now)
                return "Headers different";
            }
            if (!comparison.areEqual()) {
                return "SamComparison failed";
            }
        }

        //At this point we know that the files are not byte-wise identical, but are equal according to SamComparison and their headers are equal
        //So we iterate over reads and compare them one by one.
        try(final SamReader reader1 = getReader(sam1, validation, reference);
            final SamReader reader2 = getReader(sam2, validation, reference)) {
            final SAMRecordIterator it1 = reader1.iterator();
            final SAMRecordIterator it2 = reader2.iterator();
            while (it1.hasNext() && it2.hasNext()) {
                final SAMRecord read1 = it1.next();
                final SAMRecord read2 = it2.next();
                final String eqMessage = readsEqualAllowAddingAttributes(read1, read2);
                if (eqMessage != null){
                    return eqMessage;
                }
            }
            if (it1.hasNext() || it2.hasNext()) {
                //at least one has no more records (because the while loop is done) and at least one does have records. So we're not equal.
                return "Not the same number of reads";
            }
            return null;
        }
    }

    /**
     * Compares the reads but ignores order or attributes.
     * Also allows read2 to have a superset of attributes of read1.
     */
    private static String readsEqualAllowAddingAttributes(final SAMRecord read1, final SAMRecord read2) {
        final SAMRecordCoordinateComparator coordinateComparator = new SAMRecordCoordinateComparator();
        if (!Objects.equals(read1.getFlags(), read2.getFlags())){
            return "flags different read1:" + read1.getReadName() + " read2:" + read2.getReadName() ;
        }
        if (!Objects.equals(read1.getInferredInsertSize(), read2.getInferredInsertSize())){
            return "getInferredInsertSize different read1:" + read1.getReadName() + " read2:" + read2.getReadName() + " (" + read1.getInferredInsertSize()+ " vs " + read2.getInferredInsertSize() + ")";
        }
        if (!Objects.equals(read1.getMappingQuality(), read2.getMappingQuality())){
            return "getMappingQuality different read1:" + read1.getReadName() + " read2:" + read2.getReadName() ;
        }
        //Note the comparator does check flags, MQ and insert size but we want better messages
        if (0 != coordinateComparator.fileOrderCompare(read1, read2)){
            return "file order read1:" + read1.getContig() + ":" + read1.getAlignmentStart() + " read2:" + read2.getContig() + ":" + read2.getAlignmentStart();
        }
        if (!Objects.equals(read1.getReadName(), read2.getReadName())){
            return "name different read1:" + read1.getReadName() + " read2:" + read2.getReadName() ;
        }
        if (!Objects.equals(read1.getCigar(), read2.getCigar())){
            return "getCigar different read1:" + read1.getReadName() + " read2:" + read2.getReadName() ;
        }
        if (!Arrays.equals(read1.getReadBases(), read2.getReadBases())){
            return "getReadBases different read1:" + read1.getReadName() + " read2:" + read2.getReadName() ;
        }
        if (!Arrays.equals(read1.getBaseQualities(), read2.getBaseQualities())){
            return "getBaseQualities different read1:" + read1.getReadName() + " read2:" + read2.getReadName() ;
        }
        final List<SAMRecord.SAMTagAndValue> attr1 = read1.getAttributes();
        final List<SAMRecord.SAMTagAndValue> attr2 = read2.getAttributes();
        if (attr1.size() > attr2.size()){
            return "lost attributes read1:" + read1.getReadName() + " read2:" + read2.getReadName() + "("+ attr1.size() + " vs " + attr2.size() +")" ;
        }

        //We want to compare attributes regardless of order, so we sort first attributes by name
        final SortedMap<String, Object> values1 = new TreeMap<>();
        final SortedMap<String, Object> values2 = new TreeMap<>();

        for (SAMRecord.SAMTagAndValue samTagAndValue : attr1) {
            values1.put(samTagAndValue.tag, samTagAndValue.value);
        }
        for (SAMRecord.SAMTagAndValue samTagAndValue : attr2) {
            values2.put(samTagAndValue.tag, samTagAndValue.value);
        }

        for (int i = 0; i < attr1.size(); i++) {
            final String tag = attr1.get(i).tag;
            if (!values2.containsKey(tag)){
                return "attribute " + tag + " is missing from read2:" + read2.getReadName() ;
            }
            final Object v1 = values1.get(tag);
            final Object v2 = values2.get(tag);
            if (!Objects.equals(v1, v2)){
                return "attribute " + tag + " has different value in read1:" + read1.getReadName() + " read2:" + read2.getReadName() ;
            }
        }
        return null;
    }
}

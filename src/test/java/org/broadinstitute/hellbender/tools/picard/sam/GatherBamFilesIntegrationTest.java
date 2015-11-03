package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class GatherBamFilesIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/GatherBamFiles");
    private static final File ORIG_BAM = new File(TEST_DATA_DIR, "orig.bam");
    private static final List<File> SPLIT_BAMS_UNMAPPED_LAST = Arrays.asList(
            new File(TEST_DATA_DIR, "indchr1.bam"),
            new File(TEST_DATA_DIR, "indchr2.bam"),
            new File(TEST_DATA_DIR, "indchr3.bam"),
            new File(TEST_DATA_DIR, "indchr4.bam"),
            new File(TEST_DATA_DIR, "indchr5.bam"),
            new File(TEST_DATA_DIR, "indchr6.bam"),
            new File(TEST_DATA_DIR, "indchr7.bam"),
            new File(TEST_DATA_DIR, "indchr8.bam"),
            new File(TEST_DATA_DIR, "indUnknownChrom.bam")
    );

    private static final List<File> SPLIT_BAMS_UNMAPPED_FIRST = Arrays.asList(
            new File(TEST_DATA_DIR, "indUnknownChrom.bam"),
            new File(TEST_DATA_DIR, "indchr1.bam"),
            new File(TEST_DATA_DIR, "indchr2.bam"),
            new File(TEST_DATA_DIR, "indchr3.bam"),
            new File(TEST_DATA_DIR, "indchr4.bam"),
            new File(TEST_DATA_DIR, "indchr5.bam"),
            new File(TEST_DATA_DIR, "indchr6.bam"),
            new File(TEST_DATA_DIR, "indchr7.bam"),
            new File(TEST_DATA_DIR, "indchr8.bam")
    );


    public String getTestedClassName() {
        return GatherBamFiles.class.getSimpleName();
    }

    @Test
    public void testTheGathering() throws Exception {
        final File outputFile = BaseTest.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        final List<String> args = new ArrayList<>();
        for (final File splitBam : SPLIT_BAMS_UNMAPPED_LAST) {
            args.add("--INPUT");
            args.add(splitBam.getAbsolutePath());
        }
        args.add("--OUTPUT");
        args.add(outputFile.getAbsolutePath());
        runCommandLine(args);
        SamAssertionUtils.assertSamsEqual(ORIG_BAM, outputFile);
        Assert.assertNotNull(SamAssertionUtils.samsEqualStringent(ORIG_BAM, SPLIT_BAMS_UNMAPPED_LAST.get(0), ValidationStringency.DEFAULT_STRINGENCY, null)); // sanity check
    }

    @Test
    public void testTheGatheringUnmappedFirst() throws Exception {
        final File outputFile = BaseTest.createTempFile("gatherBamFilesTest.samFile.", BamFileIoUtils.BAM_FILE_EXTENSION);
        final List<String> args = new ArrayList<>();
        for (final File splitBam : SPLIT_BAMS_UNMAPPED_FIRST) {
            args.add("--INPUT");
            args.add(splitBam.getAbsolutePath());
        }
        args.add("--OUTPUT");
        args.add(outputFile.getAbsolutePath());
        runCommandLine(args);
        Assert.assertNotNull(SamAssertionUtils.samsEqualStringent(ORIG_BAM, outputFile, ValidationStringency.DEFAULT_STRINGENCY, null));
    }

}

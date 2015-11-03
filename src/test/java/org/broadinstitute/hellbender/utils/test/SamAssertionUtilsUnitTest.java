package org.broadinstitute.hellbender.utils.test;

import htsjdk.samtools.ValidationStringency;
import org.junit.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class SamAssertionUtilsUnitTest extends BaseTest{

    @DataProvider(name = "bamPairs")
    public Object[][] differentFilesButSameContent(){
        final String fileBam= "file1.bam";
        final String fileSam = "file1.sam";
        return new Object[][]{
                {fileBam, fileBam, true},
                {fileSam, fileSam, true},
                {fileBam, fileSam, true},
                {fileSam, fileBam, true},
                {fileBam, "file1_reorder_read_attributes.sam", true},
                {fileBam, "file1_new_read_attribute.sam", true},        //ok to add an attribute

                {fileBam, "file1_reorder_header_lines.sam", false},
                {fileBam, "file1_different_version.sam", false},
                {fileBam, "file1_missing_read_attribute.sam", false},  //not ok to lose an attribute
                {fileBam, "file1_different_bases.sam", false},
                {fileBam, "file1_different_basequals.sam", false},
                {fileBam, "file1_different_attributes.sam", false},
                {fileBam, "file1_different_mappingQ.sam", false},
                {fileBam, "file1_different_cigar.sam", false},
                {fileBam, "file1_different_position.sam", false},
        };
    }

    private static final File TEST_DATA_DIR = new File(publicTestDir, "org/broadinstitute/hellbender/utils/test/SamAssertionUtilsUnitTest");

    @Test(dataProvider = "bamPairs")
    public void testCompareStringent(final String fName1, final String fName2, boolean expectedEqual) throws Exception {
        final File f1= new File(TEST_DATA_DIR, fName1);
        final File f2= new File(TEST_DATA_DIR, fName2);
        Assert.assertEquals(expectedEqual, null == SamAssertionUtils.samsEqualStringent(f1, f2, ValidationStringency.LENIENT, null));
    }
}

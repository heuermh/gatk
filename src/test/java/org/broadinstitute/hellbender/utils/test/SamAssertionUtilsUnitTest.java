package org.broadinstitute.hellbender.utils.test;

import org.junit.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class SamAssertionUtilsUnitTest extends BaseTest{

    @DataProvider(name = "bamPairs")
    public Object[][] differentFilesButSameContent(){
        final String file1= "file1.bam";
        return new Object[][]{
                {file1, file1, true},
                {file1, "file1.sam", true},
                {"file1.sam", file1, true},
                {file1, "file1_reorder_read_attributes.sam", true},

                {file1, "file1_reorder_header_lines.sam", false},
                {file1, "file1_different_version.sam", false},
                {file1, "file1_missing_read_attribute.sam", false},
                {file1, "file1_new_read_attribute.sam", false},
                {file1, "file1_different_bases.sam", false},
                {file1, "file1_different_basequals.sam", false},
                {file1, "file1_different_attributes.sam", false},
                {file1, "file1_different_mappingQ.sam", false},
                {file1, "file1_different_cigar.sam", false},
                {file1, "file1_different_position.sam", false},
        };
    }

    private static final File TEST_DATA_DIR = new File(publicTestDir, "org/broadinstitute/hellbender/utils/test/SamAssertionUtilsUnitTest");

    @Test(dataProvider = "bamPairs")
    public void testCompareStringent(final String fName1, final String fName2, boolean expectedEqual) throws Exception {
        final File f1= new File(TEST_DATA_DIR, fName1);
        final File f2= new File(TEST_DATA_DIR, fName2);
        Assert.assertEquals(expectedEqual, SamAssertionUtils.samsEqualStringent(f1, f2));
        Assert.assertEquals(expectedEqual, SamAssertionUtils.samsEqualStringent(f2, f1));
    }
}

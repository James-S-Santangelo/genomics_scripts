# Script to create tab-delimited text file with sample IDs and paths to BAM files
# for use with Qualimap.

import os
import glob
from tqdm import tqdm

# For testing
inpath = "../../data/raw/test_validate/"
outfile = "../../resources/bamPaths_forQualimapTesting.txt"

# For real
#inpath = "../../data/raw/bam/bwa_mapping"
#outfile = "../../resources/bamPaths_forQualimap.txt"

with open(outfile, "w") as f:
    for bam in tqdm(glob.glob(inpath + "*.bam", recursive=False)):
        basename = os.path.basename(bam)
        sample = basename.split("_merged")[0]
        line = "{0}\t{1}\n".format(sample, bam)
        f.write(line)

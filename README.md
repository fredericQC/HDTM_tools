# HDTM_tools
High density transposon mutagenesis tools

This repository contains python scripts that can be used to partially process Transposon Directed Insertion Sequencing (TraDIS) data.

# sam2sites.py
This script takes as input a sam file (aligned reads) and does several things:
-Compute "true" insertion site position (the middle base pair of the 9 bp Tn5 insertion site duplication representing the insertion site)
-Compute the read count of each of those "true" sites
-Output bed and bedgraph file, stranded or not, normalized (based on the sequencing depth) or not.

# sites2genes.py
This script takes as input the output of bedtools intersect and compute gene level insertion statistics.

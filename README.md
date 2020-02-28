# HDTM_tools
High density transposon mutagenesis tools

This repository contains python scripts that can be used to partially process Transposon Directed Insertion Sequencing (TraDIS) data.

# sam2sites.py

  This script takes as input a sam file (aligned reads) and :
  - Computes "true" insertion site position (the middle base pair of the 9 bp Tn5 insertion site duplication representing the insertion site)
  - Computes the read count of each of those "true" sites
  - Outputs bed and bedgraph file, stranded or not, normalized (based on the sequencing depth) or not.

  Python related requirements:
  - python3
  - numpy
  - pandas
  - click
  
  Thanks to Audrey Bioteau (https://github.com/ABTsutenkyo) who made this script much cleaner.


# sites2genes.py

  This script takes as input the output of bedtools intersect and computes gene level insertion statistics.


  Python related requirements:
  - python3
  - numpy
  - pandas

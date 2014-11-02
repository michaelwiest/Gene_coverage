This program takes a queried gene (provided that it is in the homo_sapiens_genes.txt reference) and generates 3 plots:

1) Broken up by transcript ID, the plot displays untranslated regions (UTR), coding sequences, and coverage by Recombien (V3), Trusight, and Nextera.

2) Shows the same information as 1 but not broken up by transcript ID. The different transcript IDs are overlaid and semi-transparent.

3) Similar to 2 except the transcript IDs are staggered vertically for increased readability.


Also a text file is generated that displays the different coverage regions for the reference file and for the different sequencing technologies.
If a sequencing technology displays all zeros then it has no coverage for the given gene. 
_____________________________________

Dependancies: Python: numpy
	      R: plyr, dplyr, gpplot2, data.table


_____________________________________

usage:

./GeneComparator.sh <somegene>

eg: ./GeneComparator.sh HBB


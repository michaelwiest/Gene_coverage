


This program takes a queried gene (provided that it is in resources/homo_sapiens_genes.txt reference) and generates 3 plots:

1) Broken up by transcript ID, the plot displays untranslated regions (UTR), coding sequences, and coverage by Recombine (V3), Trusight, and Nextera.

2) Shows the same information as 1 but not broken up by transcript ID. The different transcript IDs are overlaid and semi-transparent.

3) Similar to 2 except the transcript IDs are staggered vertically for increased readability.


Also a text file is generated that displays the different coverage regions for the reference file and for the different sequencing technologies.
If a sequencing technology displays all zeros then it has no coverage for the given gene. 
_____________________________________

Dependancies: 
	      Python: numpy
	      R: gpplot2, data.table


_____________________________________

USAGE:

./GeneCoverage.sh <somegene>

eg: ./GeneCoverage.sh HBB

HBB is a good gene to query because it is very easy to see the various coverages. Very large genes can be hard to view and use of the *.txt file is often preferable.





NOTE: The reference file for this project is HG37, not HG38.


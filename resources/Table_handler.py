import sys, os, numpy, re

rawdata=sys.argv[1]
outfolder="resources"



#READ THE SPECIFIED COLUMNS FROM THE RAW TEXT FILE
#THEY ARE THE CHROMOSOME, FEATURE, START, AND END
table=numpy.loadtxt(rawdata, unpack=False, dtype='str', usecols = (0,2,3,4))

#GENERATE THE CORRECT PATH
tableout=os.path.join(outfolder, "TEMP.gene.table.txt")

numpy.savetxt(tableout, table, fmt='%s')

#THIS SECTION IS FOR EXTRACTING TRANSCRIPT IDS
table2=rawdata

transcript_id=[]
with open(table2) as input:
	for lines in input:
		lines=str(lines)
		#FIND THE TRANSCRIPT ID
		tid=re.search('transcript_id "(.+?)";', lines)
		if tid:
			tid=tid.group(1)
			transcript_id.append(tid)
		else:
			transcript_id.append("NA")

tidout = os.path.join(outfolder, "TEMP.gene.tid.txt")

#OUTPUT THE TRANSCRIPT IDS AS AN ARRAY TO MAKE IT READABLE TO R
numpy.savetxt(tidout, transcript_id, fmt='%s')


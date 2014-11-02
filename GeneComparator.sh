#!/bin/bash

function execute {

#CHECK THAT THE QUERIED GENE IS IN THE HOMO SAPIENS FILE
#THIS DOES NOT CHECK THAT IT IS IN ANY OF THE OTHER CHIP/AMPLICON FILES
if grep -q $1 "resources/homo_sapiens_genes.txt"; then

mkdir -p output
echo "
Extracting lines from HomoSapiens file...
"

grep -i  "\<$1\>" resources/Homo_sapiens.GRCh37.75.gtf > resources/TEMP.gene.raw.txt

echo "
Processing data in Python...
"

Python resources/Table_handler.py resources/TEMP.gene.raw.txt

echo "
Handling data and printing plots with R...
"

Rscript --vanilla  resources/Gene_plotter.5.R resources/TEMP.gene.table.txt resources/TEMP.gene.tid.txt $1


rm resources/TEMP.*.txt


else
echo "

Please query a valid gene. See resources/homo_sapiens_genes.txt for list of genes.

"


fi

}
#CHECK THAT THE HOMOSAPIENS FILE IS PRESENT.

if [ -e "resources/Homo_sapiens.GRCh37.75.gtf" ]; then
echo "Homo_sapiens file present, beginning..."
execute $1

else
while true; do
read -p "Do you want to download the Homo_sapiens resource GTF file? This is large (>100MB), but is necessary to run. (y/n)  " yn

#GET THE HOMO SAPIENS FILE
case $yn in 
[Yy]* ) curl -O ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz 
echo "unzipping file..."
gunzip "Homo_sapiens.GRCh37.75.gtf.gz" ; 
mv "Homo_sapiens.GRCh37.75.gtf" "resources" ; 
execute $1 ; exit ;; 

[Nn]* ) echo "quitting..." ; exit 1;;
* ) echo "Please answer yes or no." ;;

esac
done



fi



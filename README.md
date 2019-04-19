# BifIDaC

BifIDaC (Bifidobacterial Its Database Creator) is a python tool to create an ITS amplicon database specific of genus
Bifidobacterium

BifIDac requires Python 3.x, Biopython library and a working version of Clustalw installed. 

Folder comes with two files, called fw_primer and rev_primer.
By default they contain sequence of  Probio-bif_Uni and Probio-bif_Rev , designed and tested in 
"Evaluation of bifidobacterial community composition in the human gut by means of a targeted 
amplicon sequencing (ITS) protocol" by Milani C. et al., 2014. 
https://doi.org/10.1111/1574-6941.12410

Header of sequences is not relevant. If namefile is changed, adjustment needs to be done
into code.
Every retrieved sequence was written in direction 5' -> 3'
Header contains the associated taxonomy.


Enjoy.


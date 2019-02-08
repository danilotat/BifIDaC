#!/usr/bin/env python3

"""

	BIF_its is a python script written to make easier the making of .fasta database composed by
	every amplicon generated from a couple of defined primers Probio-bif_Uni (5′-CTKTTGGGYYCCCKGRYYG-3′)
	Probio-bif_Rev (5′-CGCGTCCACTMTCCAGTTCTC-3′) inside the genus Bifidobacterium.
	
	This POC code is just an experiment made for lulz to practice with biopython modules and some 
	filetype used in bioinformatics, like genbank, fasta, and .aln file produced by Clustalw. This one
	could be easly extended to other genus, other primers: by far, this will fit properly with any kind of 
	intergenic extraction or amplicon selection.
	
	Code was tested on:	 Windows10, x64 with Python 3.x
				 Ubuntu 18.04 LTS, x64 with Python 3.x 
	
	I don't give a shit about lacking pythonic syntax.
	
	Made for fun and practice by @danilotat
	
"""
	
	
import os
import shutil
import glob
import platform
import sys
import importlib.util
import time
from Bio import SeqIO
from Bio import Entrez
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

# check presence of biopython

package = "Bio"
spec = importlib.util.find_spec(package)
print(spec)
if spec is None:
    print(package+" is not installed.")
    print("Please install using")
    print("$ pip -install biopython")
    sys.exit()

# check right existance of primers file

cwd = os.getcwd()
primer1 = cwd+'/rev_primer.fasta'
primer2 = cwd+'/fw_primer.fasta'
if not os.path.isfile(primer1):
    print('Reverse primer seems to miss.')
    print('Quitting..')
    sys.exit()
if not os.path.isfile(primer2):
    print('Forward primer seems to miss.')
    print('Quitting..')
    sys.exit()


# path of Clustalw + check existance


clustalw = 0
platform = platform.system()
if platform == "Linux":
    clustalw = "/usr/bin/clustalw"
    if not os.path.isfile(clustalw):
        print("Clustalw is not installed.")
        print("Please install using $ sudo apt-get install clustalw ")
        sys.exit()
elif platform == "Windows":
    clustalw = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    if not os.path.isfile(clustalw):
        print("Clustalw is not installed.")
        print("Please download from the official website")
        sys.exit()

#folder setup
temp_folder = str(cwd)+'/temp'
if not os.path.exists(temp_folder):
    print('Creating new temp folder..')
    os.makedirs(temp_folder)
else:
    print('Temp folder already exists. Removing..')
    shutil.rmtree(temp_folder)
    os.makedirs(temp_folder)

gbk_folder = str(cwd)+'/genomes'
if not os.path.exists(gbk_folder):
    print('')
    print("Creating a new folder..   ")
    print('')
    os.makedirs(gbk_folder)
else:
    print('Folder already exists. Removing older files..   ')
    shutil.rmtree(gbk_folder)
    os.makedirs(gbk_folder)

#   genome downloading
#   time sleep is used to avoid server timeout while querying

Entrez.email = ""	#insert email
Entrez.api_key = ""	#insert your api key
search_term = "Bifidobacterium[organism] AND complete+genome[title]"
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=500)
genome_ids = Entrez.read(handle)['IdList']
print(genome_ids)
for genome_id in genome_ids:
    record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gbwithparts", retmode="text")
    filename = 'GenBank_Record_{}.gbk'.format(genome_id)
    print('Writing:{}'.format(filename))
    with open(filename, "w") as f:
        f.write(record.read())
    time.sleep(0.5)


# moving genomes into directory

genomes = glob.glob(cwd+'/*.gbk')
for file in genomes:
    shutil.move(file,gbk_folder)

print("Getting sequences from selected genbank file..  ")
start = 0
end = 0
strand = 0

output_file = temp_folder+"/ITS_sequence.fasta"
intergenic_sequence = []
for f in glob.glob(cwd+'/genomes/*.gbk'):
    gbank = SeqIO.parse(open(f,'r'), "genbank")
    with open(output_file, "w") as nfh:
        for genome in gbank:
            # trait for ITS located onto reverse strand
            rnas_rev = [feat for feat  in genome.features if feat.type == 'rRNA' and feat.location.strand == -1]
            s_rev = [ r for r in rnas_rev if '16S' in r.qualifiers['product'][0]]
            lo_rev = [ r for r in rnas_rev if '23S' in r.qualifiers['product'][0]]
            # trait for ITS located onto forward strand
            rnas_fw = [ feat for feat  in genome.features if feat.type == 'rRNA' and feat.location.strand == 1]
            s_fw = [ r for r in rnas_fw if '16S' in r.qualifiers['product'][0]]
            lo_fw = [ r for r in rnas_fw if '23S' in r.qualifiers['product'][0]]
                # i select 16S
                # j select 23S
            for i,j in zip(s_rev,lo_rev):
                start = j.location.nofuzzy_end
                end = i.location.nofuzzy_start
                ig_seq = genome.seq[start:end]
                ig_seq_rev_compl = ig_seq.reverse_complement()
                intergenic_sequence.append(">%s|ITS_full_sequence|\n%s" % (genome.description, ig_seq_rev_compl))
            for i,j in zip(s_fw,lo_fw):
                start = i.location.nofuzzy_end
                end = j.location.nofuzzy_start
                intergenic_sequence.append(">%s|ITS_full_sequence|\n%s" % (genome.description, genome.seq[start:end]))
        for item in intergenic_sequence:
            nfh.write("%s\n" % item)

# delete empty sequence from file
print("Removing empty sequences.. ")
no_empty_seq = temp_folder+'/no_empty_seq.fasta'
seqs = []
with open(no_empty_seq, "w") as nk:
    for record in SeqIO.parse(output_file, "fasta"):
        if len(record.seq):
            seqs.append(record)
SeqIO.write(seqs, no_empty_seq, "fasta")

# primer joining

print("Joining files..")

out_file_1 = temp_folder+"/first_assembled.fasta"
file = [cwd+"/fw_primer.fasta", no_empty_seq]
sep = ','
sep2 = ' ='
with open(out_file_1, "w") as o_file_1:
    for fname in file:
        with open(fname) as infile:
            for line in infile:
                if line.startswith(">"):
                    short_line = line.split(sep,1)[0]
                    short_line1 = short_line.split(sep2, 1)[0]
                    new_line = short_line1.replace(' ', '_')
                    new_line_2 = new_line.replace('Bifidobacterium_', 'B.')
                    new_line_3 = new_line_2.replace('strain_', '')
                    new_line_4 = new_line_3.replace('subsp._', 's.')
                    new_line_5 = new_line_4.replace('_complete_genome', '')
                    new_line_6 = new_line_5.replace('_chromosome', '')
                    new_line_7 = new_line_6.replace(',', '')
                    o_file_1.write(new_line_7)
                    o_file_1.write("\n")
                else:
                    o_file_1.write(line)

# Duplicate remover

print("I'm removing duplicates..")
idlist = set()
seqiter = SeqIO.parse(out_file_1, 'fasta')
with open(temp_folder + '/ITS_noduplicate.fasta', 'w') as nhl:
    for record in seqiter:
        if record.id not in idlist:
            idlist.add(record.description)
            SeqIO.write(record,nhl,'fasta')
print("Duplicates are removed")

# alignment

clustalw_in_1 = temp_folder+"/ITS_noduplicate.fasta"
clustaw_out_1 = temp_folder+"/1st_aln"
clustal_cline = ClustalwCommandline(clustalw, infile=clustalw_in_1, outfile=clustaw_out_1)
assert os.path.isfile(clustalw), "Clustal W executable missing"
stdout, stderr = clustal_cline()

# trim on alignment based on primer fw align profile
print('Alignment..')
output = temp_folder+"/experiment.fasta"
with open(output, "w") as gls:
    aln = AlignIO.read(clustaw_out_1, "clustal")
    for col in range(aln.get_alignment_length()):
        res = list(aln[:,col])
        if not '-' in res:
            position = col
            print('First full column is {}'.format(col))
            first_alignment = (aln[:,position:])
            SeqIO.write(first_alignment, gls, "fasta")
            break

# remove gaps and primer sequence

output_file_2 = temp_folder+"/aln_second_edit.fasta"
with open(output, "r") as nfh:
    with open(output_file_2, "w") as nk:
        for line in nfh:
            if not '>' in line:
                no_gaps_seq = line.replace('-', '')
                nk.write(no_gaps_seq)
            else:
                nk.write(line)


with open(temp_folder+"/aln_third_edit.fasta", "w") as ekr:
    edit = SeqIO.parse(output_file_2, "fasta")
    for record in edit:
        if len(record.seq) > 50:
            SeqIO.write(record, ekr, "fasta")


# adding forward primer to aligned file and do same things as above

output_file_2 = temp_folder+"/aln_forth_edit.fasta"
with open(output_file_2, "w") as ekw:
    file=[cwd+"/rev_primer.fasta", temp_folder+"/aln_third_edit.fasta"]
    for fname in file:
        with open(fname) as kil:
            for line in kil:
                ekw.write(line)

# 2nd alignment

print("2nd alignment..")
clustaw_out_2 = temp_folder+"/2nd_aln"
clustal_cline=ClustalwCommandline(clustalw, infile=output_file_2, outfile=clustaw_out_2)
assert os.path.isfile(clustalw), "Clustal W executable missing"
stdout, stderr = clustal_cline()

# last alignment edit

output=temp_folder+"/2nd_aln_1st.fasta"
with open(output, "w") as gls:
    aln = AlignIO.read(clustaw_out_2, "clustal")
    for col in range(aln.get_alignment_length()):
        res = list(aln[:, col])
        if not '-' in res:
            position = col + 21     #lenght of other primer
            print('First full column is {}'.format(col))
            first_alignment = (aln[:, :position ])
            SeqIO.write(first_alignment, gls, "fasta")
            break

with open(temp_folder+ '/ITS_seq_no_gaps.fasta', "w") as nfh:
    with open(output, "r") as nkl:
        for line in nkl:
            new_line = line.replace('-', '')
            nfh.write(new_line)

with open(cwd + '/ITS_final.fasta', 'w') as nfh:
    edit = SeqIO.parse(temp_folder + '/ITS_seq_no_gaps.fasta', 'fasta')
    for record in edit:
        if len(record.seq) > 50:
            SeqIO.write(record, nfh, "fasta")

if os.path.exists(cwd+'/temp'):
    shutil.rmtree(cwd+'/temp')
if os.path.exists(cwd+'/genomes'):
    shutil.rmtree(cwd+'/genomes')

print('')
print('Work is done')
print("You'll find the ITS database file under the current directory")
print("Report any error to @danilotat")

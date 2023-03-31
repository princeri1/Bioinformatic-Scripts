from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Phylo
import sys
import os
import subprocess
import re

# this section takes care of reading in data from user
usage = "Usage: Project.py <DNA Sequence fasta file>"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

# read in the fasta file name and build an Ascii tree
fasta = sys.argv[1]
outfile = fasta[0:len(fasta)-6:1]
os.system("clustalw2 -infile=./"+ fasta + " -type=DNA -outfile=./"+ outfile + ".aln")
tree = Phylo.read(outfile+".dnd","newick")
Phylo.draw_ascii(tree)



# parse the fasta file and seperate into sequences
data = SeqIO.parse(fasta,"fasta")
sequences = {}
keys = []
for line in data :
	keys.append(line.id)
	sequences.update({line.id : line.seq})
#Creates a folder to store individual RNA Sequences
os.system("mkdir RNAfiles")

#Creates Fasta files for each Sequence and Runs RNAfold on each of them
foldOut = {}
for i in keys :
	RNA = sequences[i].transcribe()
	#Create File for Each RNA Sequence
	f = open("RNAfiles/"+ i + ".fasta","w")
	f.write(">" + str(i) + "\n")
	f.write(str(RNA))
	f.close()
	#Run RNA fold on Each File Created
	filepath = "./RNAfiles/"+ i + ".fasta"
	print(filepath)
	pipeIn  = subprocess.Popen(('cat', filepath), stdout=subprocess.PIPE)
	output = subprocess.check_output(("RNAfold"), stdin=pipeIn.stdout)
	data = output.decode()
	data = data.replace('\n', ' ')
	foldOut.update({i:data})


#This prints the %of nucleotides that are folded into hairpins
#data Stored incase needed for future applications
data = []
print("Sequence Name, Percent of folded nucleotides, Energies")
for i in keys :
	name, sequence, brackets, energy = foldOut[i].split()
	percent_matched = (2 * brackets.count("(")) / (len(sequence))
	tupe = (i,percent_matched, energy)
	print(tupe)
	data.append(tupe)

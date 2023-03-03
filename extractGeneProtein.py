#Riley Prince BB452

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import re

# this section takes care of reading in data from user
usage = "Usage: " + sys.argv[0] + " <genome fasta> <gff file> <transcript ID>"
if len(sys.argv) != 4:
    print usage
    sys.exit()

# read the input files/args.
genomeFasta = sys.argv[1]
gffFile = sys.argv[2]
name = sys.argv[3]

# read the fasta file into a dictionary.
genome = {}
sequences = SeqIO.parse(genomeFasta,"fasta")
for record in sequences:
    genome[record.id] = record.seq

# initialize some variables
revcomp = False
coords = []

# read through the gffFile          
GFF = open(gffFile,'r')
for line in GFF:
    if "#" not in line:
        chrom,source,seqtype,start,stop,score,strand,frame,attributes = line.strip().split("\t")
        if "CDS" in seqtype:
            if name in attributes:
                coords.append((chrom,int(start)-1,int(stop)))
                if strand == "-":
                    revcomp = True

# reverse the order of the exons if on "-" strand
coords.sort(reverse = revcomp)

# collect the coding exons of the transcript 
ORF = Seq('')
for chrom,start,stop in coords:
    CDS = genome[chrom][start:stop]
    if revcomp:
        CDS = genome[chrom][start:stop].reverse_complement()
    # concatenate the coding sequence 
    ORF += CDS

# transcribe,translate, and print 
RNA = ORF.transcribe()
Protein = RNA.translate()
print Protein

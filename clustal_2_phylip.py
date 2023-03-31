import sys
from Bio import AlignIO

# this section takes care of reading in data from user
usage = "Usage: clustal_2_phylip.py <.aln file>"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

# read in the fastq file name
aln_file = sys.argv[1]
phy_file = aln_file[:len(aln_file)-3:1] + ".php"

alignment = AlignIO.parse(aln_file,"clustal")
AlignIO.write(alignment,open("18S_rRNAs.phy","w"),"phylip")



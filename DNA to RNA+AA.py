#Riley Prince BB452
from Bio import SeqIO
from Bio import Seq
from Bio import SeqUtils
import sys

usage = "Usage: " + sys.argv[0] + " <Sequence Fasta> "
if len(sys.argv) != 2:
    print usage
    sys.exit()

fastaFile = sys.argv[1]

sequences = SeqIO.parse(fastaFile,"fasta")
instances = []

for data in sequences :
	instance = data.seq
	instances.append(instance)

DNA = Seq.Seq(instances[0])
mRNA = DNA.transcribe()
AA = mRNA.translate()

print("mRNA: ", mRNA)
print("\n", "Amino Acid Chain: ", AA)


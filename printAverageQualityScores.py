#Riley Prince bb452
import sys
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# this section takes care of reading in data from user
usage = "Usage: printAverageQualityScores.py <fastq file>"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

# read in the fastq file name
fastq = sys.argv[1]
# parse the fastq file
data = SeqIO.parse(fastq,"fastq")

# initialize a list of probabilities 
sum_p = [0] * 200
sum_q = [0] * 200
N = [0] * 200
for record in data:
    for i,Q in enumerate(record.letter_annotations["phred_quality"]):
        # convert the PHRED score to a probability
        p_err = 10**(-float(Q)/10)
        # append this specific probabity to array for this sequence
        sum_q[i] += Q
        sum_p[i] += p_err
        N[i] += 1
# now for plotting. Initialize x and y arrays to plot:
x = []
y = []

# add the average probabilties to the y values:
for i in range(len(N)):
    if N[i] > 0:
        qAvg = sum_q[i]/N[i]
        x.append(i)
        y.append(qAvg)
        
# plot the x and y values:
plt.plot(x,y, label="Prince")
plt.title('Avg Quality Score vs Position')
plt.xlabel('position (nt)')
plt.ylabel('Averge Phred Quality Score')
plt.legend(loc="upper left")
plt.savefig('quality_avg.png')

from Bio import Seq
import random

length = 200

random.seed()
i = 0
sequence = ""
#Build Palindrome Sequence

while i < length / 2 :
	randomness = random.randint(1,4)
	if randomness % 4 == 0 :
		sequence += "A"
	elif randomness % 3 == 0 :
		sequence += "C"
	elif randomness % 1 == 0 :
		sequence += "G"
	else :
		sequence += "T"
	i = i + 1

DNA = Seq.Seq(sequence)
palindrome = (DNA + DNA.reverse_complement())
if palindrome == palindrome.reverse_complement() :
	isPalindrome = "Yes!"
else :
	isPalindrome = "No :("
print(len(palindrome))
print("Your randomly generated palindrome is:\n", palindrome, "\n")
print("Is this a Palindrome?  " + isPalindrome)


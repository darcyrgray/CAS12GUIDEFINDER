# New program to take Mec guides and check against multiple SA genomes to check in all
# Input meconly_SCCmec_guide_RNAs.txt - output file of CAS12guidefinder_sa.py

from Bio.Seq import Seq
from Bio import SeqIO

print("Checking lines against MRSA genomes")

guide_dict = {};    #Dictionary used to increment if guide present in local genome    

input =  open("meconly_SCCmec_guide_RNAs.txt", "r"); #Opens guides - output at TXT file from CAS12guidefinder.py

# parse input file into above dictonaries 
for line in input:
	line=line[0:25] # 25 bp long crRNA guides
	guide_dict[line]=0;
# Function used to count occurances in a sequence string of each guide - sequence record must be changed into string for this
def count_sequence(crRNA, seq):
	return seq.count(crRNA)


#Now LVH genomes of SA genomes as FASTA files
genomes = SeqIO.parse("new_concatenated_contigs.fasta","fasta");
countg = 0;

listgenomes=[];
# iterate over genomes
for seq_record in genomes:
	countg+=1; #incremental counter to count how many genomes in fasta record checked
	listgenomes.append(seq_record.description);
	checkseq = seq_record.seq + seq_record.seq.reverse_complement();


	for guide in guide_dict:
		if guide in checkseq:
			guide_dict[guide] +=1;



#Print output to TXT file for review

out = open("MEC_GUIDES-checked.txt","wt");

for guide in guide_dict:
	percent = (guide_dict[guide]/countg)*100
	out.write(f"++++ {guide} ++++\n");
	out.write(f"Present in {percent} % of checked genomes\n");
	out.write("++++ ++++ ++++ ++++ +++++ +++++ ++++ ++++ \n\n")


out.write(f"Checked {countg} MRSA genomes\n");
print("Finished");




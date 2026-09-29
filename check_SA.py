# New program to take SA guides and check against multiple SA genomes to check in all
# Input CAS12_guide_RNAs.txt - output file of CAS12guidefinder_sa.py

from Bio.Seq import Seq
from Bio import SeqIO

print("Checking lines against SA genomes")

guide_dict = {};    #Dictionary used to increment if guide present in reference genome    
guide_dict_lv ={};  # Dictionary used to increment if guide present in local genome
guide_counts = {} #Number of copies of guide sequence in the original genome used
list_counts=[];  # will be used to count number of each guide in a genome
cumulative_count={} # Used to increment total count of copies in each reference genome - used to calculate avg copy #
cumulative_count_lv={}; #Same as line above but for local genomes 

input =  open("recheckedSAguides01AUG2026.txt", "r"); #Opens guides - output at TXT file from CAS12guidefinder.py

# parse input file into above dictonaries 
for line in input:
	linecount=line[26:30] 
	linecount = linecount.replace("\n","")
	line=line[0:25] # 25 bp long crRNA guides
	guide_dict[line]=0;
	guide_dict_lv[line]=0;
	guide_counts[line] = linecount;
	cumulative_count[line]=0;
	cumulative_count_lv[line]=0;

# use to count occurances in a sequence string of each guide - sequence record must be changed into string for this
def count_sequence(crRNA, seq):
	return seq.count(crRNA)

#Now load heaps of SA genomes as FASTA files
genomes = SeqIO.parse("SAgenomes.fasta","fasta");
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

		checkseqstr = str(checkseq);
		guidestr = str(guide);
		count_guideinstr = count_sequence(guidestr, checkseqstr);
		cumulative_count[guide]+=count_guideinstr;

# Load LVH SA Genomes as FASTA File

lvgenomes = SeqIO.parse("new_concatenated_contigs.fasta", "fasta");
countg_l = 0;

# iterate over these genomes
for record in lvgenomes:
	countg_l +=1;
	checkseq_lv = record.seq + record.seq.reverse_complement();
	for guide in guide_dict_lv:

		if  guide in checkseq_lv:
			guide_dict_lv[guide] +=1;
		checkseqstr_lv = str(checkseq);
		guidestr_lv = str(guide);
		count_guideinstr_lv = count_sequence(guidestr_lv, checkseqstr_lv);
		cumulative_count_lv[guide]+=count_guideinstr_lv;



#Print output to TXT file for review

out = open("NEW_SA_GUIDES-checked.txt","wt");

for guide in guide_dict:
	percent = (guide_dict[guide]/countg)*100
	percent_lv = (guide_dict_lv[guide]/countg_l)*100;
	avgcount = (cumulative_count[guide]/countg);
	avgcount_lv = (cumulative_count_lv[guide]/countg_l);
	out.write(f"++++ {guide} ++++\n");
	out.write(f"Present in {percent} % of checked genomes\n");
	out.write(f"Present in {guide_counts[guide]} copies of original genome\n")
	out.write(f"Present at an average copy number of {avgcount} in all included reference genomes\n");
	out.write("++++ ++++     +++++      ++++ ++++ \n")
	out.write(f"Present in {percent_lv} % of LVH MRSA genomes\n")
	out.write(f"Present at an average copy number of {avgcount_lv} in all LVH MRSA genomes\n");
	out.write("++++ ++++ ++++ ++++ +++++ +++++ ++++ ++++ \n\n")


out.write(f"Checked {countg} refseq genomes\n");

out.write("Refseq genome records checked:")

for i in listgenomes:
	out.write(f"{i}\n");

print(cumulative_count)
print(cumulative_count_lv)
print("Finished");




#Program version 2.0 #
#This program is to find CAS12 guides for a specific organism genome and then exclude primers which match (+/-1 bp wobble) those in other submitted genomes
#Finds 21 base pair sequences following PAM sequence from a given genome
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import combinations, product

#define PAM sequence to use
PAM = ["TTTC", "TTTG", "TTTA"]
BASES = "ATCG"
#####################################################################
def generate_mismatch_variants(seq, positions_range, max_mismatches):
	seq=list(seq)
	candidate_positions= list(positions_range)
	for k in range (1, max_mismatches+1):
		for pos_combo in combinations(candidate_positions, k):
			sub_options = []
			for pos in pos_combo:
				original=seq[pos]
				sub_options.append([b for b in BASES if b != original])
			for subs in product(*sub_options):
				variant = seq[:]
			for pos,base in zip(pos_combo,subs):
				variant[pos]=base
			yield "".join(variant)
####################################################################
#Identify target sequence to open
input_genome = input("Please enter full path to target genome:")

target = open(input_genome,"r")

#process genomes
targetp = target.read().replace("\n", "")
targetpf = Seq(targetp)
targetprc = Seq(targetpf).reverse_complement()

#find primer run - TTV PAM at start (C/G/A) then 18-22 further BP (aiming 21)
a = 0
posprim = ""
primer = ""

#Look for multicopy possible gRNA with dictonary incremented with occurance
dict={}

#Search forward strand #
while a < len(targetp):

	##loading progress bar
	if a%1000==0:
		progress = "#"*int(100*a/len(targetp))
		remaining = " "*int(100-(100*a/len(targetp)))
		print(f"FORWARD LOOP[{progress}+{remaining}]")

	#####
	if targetpf[a:a+4] in PAM:
		posprim=targetpf[a:a+25]
		if posprim in dict:
			dict[posprim]+=1
			a+=25
		else:
			dict[posprim]=1
			primer = "1"
			a+=25
	else:
		a+=1
		posprim= ""

#Search reverse compliment strand # 
a = 0

while a < len(targetp):

	#loading progress bar
	if a%1000==0:
		progress = "#"*int(100*a/len(targetp))
		remaining = " "*int(100-(100*a/len(targetp)))
		print(f"REVERSE LOOP[{progress}+{remaining}]")

	####
	if targetprc[a:a+4] in PAM:
		posprim=targetprc[a:a+25]
		if posprim in dict:
			dict[posprim]+=1
			a+=25
		else:
			dict[posprim]=1
			primer = "1"
			a+=25
	else:
		a+=1
		posprim= ""

#IF no primers found
if primer == "":
	print("No definite primer found")


##### Find highest copy number primers ####
#find Loopmax #  most frequent primers
maxdict = {}
n=0
loopmax = 100
while n<loopmax:
	maxcp = max(dict, key=dict.get)
	maxdict[maxcp] = dict[maxcp]
	dict[maxcp]=0

	n+=1


#Program pauses here significantly loading human genome thus the print statement
print("++++++++  NOW LOADING NON-TARGET GENOMES TO EXCLUDE OFF TARGET gRNA BINDING - KINDLY WAIT +++++++++")
exclude_file = input("Please enter full path to list of excluded PAM sequences - use exclusion PAM generator with reference genomes to be excluded: ")
exclude_file_open = open(exclude_file,"r")

exclude_list = set(exclude_file_open.read().splitlines())

#### LOOP TO EXCLUDE ###
#loop counter lOoP
lOoP = 0
MAX_MISMATCHES = 3 


for item in maxdict:
	item_str = str(item)
	if item_str in exclude_list:
		maxdict[item]=0
		continue

	#loading progress bar
	lOoP+=1
	progress = "#"*int(100*lOoP/loopmax)
	remaining = " "*int(100-(100*lOoP/loopmax))
	print(f"Exclusion loops[{progress}+{remaining}]")


	excluded = False
	for variant in generate_mismatch_variants(item_str, range(4,25), MAX_MISMATCHES):
		if variant in exclude_list:
			excluded = True
			break
		if excluded:
			break
	if excluded:
		maxdict[item]=0


#### Write guides into .txt file ###
outputfile = input("Please enter name of output file (txt):")
out = open(outputfile, "wt")

for item in maxdict:
	out.write(f"{item} {maxdict[item]}\n")
out.close()

#This program is to find CAS12 guides for a specific organism genome and then exclude primers which match (+/-1 bp wobble) those in other submitted genomes
#Finds 21 base pair sequences following PAM sequence from a given genome
from Bio.Seq import Seq
from Bio import SeqIO


#define PAM sequence to use
PAM = ["TTTC", "TTTG", "TTTA"];

#Identify target sequence to open
input_genome = input("Please enter full path to target genome:");

target = open(input_genome,"r"); 

#process genomes
targetp = target.read().replace("\n", "");
targetpf = Seq(targetp); 
targetprc = Seq(targetpf).reverse_complement();

#find primer run - TTV PAM at start (C/G/A) then 18-22 further BP (aiming 21)
a = 0;
posprim = "";
primer = "";

#Look for multicopy possible gRNA with dictonary incremented with occurance
dict={};

#Search forward strand #
while a < len(targetp):
 
	##loading progress bar

	progress = "#"*int(100*a/len(targetp))
	remaining = " "*int(100-(100*a/len(targetp)))
	print(f"FORWARD LOOP[{progress}+{remaining}]");

	#####
	if targetpf[a:a+4] in PAM:
		posprim=targetpf[a:a+25];
		if posprim in dict:
			dict[posprim]+=1;
			a+=25;
		else:
			dict[posprim]=1;
			primer = "1";
			a+=25;
	else:
		a+=1;
		posprim= "";

#Search reverse compliment strand # 
a = 0;

while a < len(targetp):
	
	#loading progress bar

	progress = "#"*int(100*a/len(targetp))
	remaining = " "*int(100-(100*a/len(targetp)))
	print(f"REVERSE LOOP[{progress}+{remaining}]");

	####
	if targetprc[a:a+4] in PAM:
		posprim=targetprc[a:a+25];
		if posprim in dict:
			dict[posprim]+=1;
			a+=25;
		else:
			dict[posprim]=1;
			primer = "1";
			a+=25;
	else:
		a+=1;
		posprim= "";

#IF no primers found
if primer == "":
	print("No definite primer found");



##### Find highest copy number primers ####
#find Loopmax #  most frequent primers
maxdict = {};
n=0;
loopmax = 20;
while n<loopmax:
	maxcp = max(dict, key=dict.get);
	maxdict[maxcp] = dict[maxcp];
	dict[maxcp]=0;

	n+=1;

#Program pauses here significantly loading human genome thus the print statement
print("++++++++  NOW LOADING NON-TARGET GENOMES TO EXCLUDE OFF TARGET gRNA BINDING - KINDLY WAIT +++++++++")
#### EXCLUDE PRESENCE of these highest copy number crRNA IN NOT TARGET ORGS ####
#Can add genomes to exclude:
exclude = "";
exclude = Seq(exclude);

genomes_exclude = [];
while True:
	Add_exclude = input("Please enter full path to genomes which should have presence of guides excluded - must be fasta formatted) -- Type 'exit' to exit when finished:")
	if Add_exclude.upper() == "EXIT": 
		break
	else:
		genomes_exclude.append(Add_exclude)
nni=0;
nameindex = {};

for item in genomes_exclude:
	name = "exseq" +f"{nni}"
	nameindex[name] = SeqIO.parse(item, 'fasta');

exclude = Seq("");
############ LOOP TO EXCLUDE ###################
for item in nameindex:
	for seq_record in nameindex[item]:
		exclude += (seq_record.seq);

excluderv = exclude.reverse_complement();

#loop counter lOoP
lOoP = 0;


for item in maxdict:

	#loading progress bar
	lOoP+=1;
	progress = "#"*int(100*lOoP/loopmax);
	remaining = " "*int(100-(100*lOoP/loopmax));
	print(f"Exclusion loops[{progress}+{remaining}]");
	####
	
#checking guides for wobble specificity (1bp) at any position afer PAM to reference CONS or Human genome (hs1.fa)
#May not be needed given checking guides with blast identifies any wobble specificity issue 
#Additionally this SIGNIFICANTLY increases compute time 

	wobble = ["A", "T", "C", "G"];
	for base in wobble:
		i=3
		#print(base); #old print statement used to monitor progress as the Wobble checking is SLOW
		while i<25:
			item_temp = item[:i]+base+item[(i+1):]; 
			i+=1;
			#print(i); # same as above - old print statment 
			if item_temp in exclude or item in excluderv:

				dict[item] = 0;



#### Write guides into .txt file ###
outputfile = input("Please enter name of output file (txt):")
out = open(outputfile, "wt");

for item in maxdict:
	out.write(f"{item} {maxdict[item]}\n");

out.close();

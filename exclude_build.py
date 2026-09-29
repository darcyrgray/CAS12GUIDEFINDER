from Bio.Seq import Seq
from Bio import SeqIO
import re 

PAM = ["TTTC", "TTTG", "TTTA"]
#Build a file - comma separated list of PAM lead 25 bp sequeuences in genomes to exclude:
#### EXCLUDE PRESENCE of these highest copy number crRNA IN NOT TARGET ORGS ####
#Can add genomes to exclude:
exclude = ""
exclude = Seq(exclude)

genomes_exclude = []
while True:
	Add_exclude = input("Please enter full path to genomes which should have presence of guides excluded - must be fasta formatted) -- Type 'exit' to exit when finished:")
	if Add_exclude.upper() == "EXIT": 
		break
	else:
		genomes_exclude.append(Add_exclude)
nni=0
nameindex = {}

for item in genomes_exclude:
	name = "exseq" +f"{nni}"
	nameindex[name] = SeqIO.parse(item, 'fasta')
	nni+=1;

exclude = Seq("")
seq_parts =[]
for item in nameindex:
	for seq_record in nameindex[item]:
		seq_parts.append(str(seq_record.seq))

exclude = Seq("".join(seq_parts))

excluderv = exclude.reverse_complement()

### Make list of guide sequences in exlcuded genomes ###
exclude_list = set()

pam_regex = re.compile(r'(?=(TTT[ACG][ACGT]{21}))')

def find_guides(seq_str):
	for m in pam_regex.finditer(seq_str):
		yield m.group(1)

for genome_seq in [str(exclude), str(excluderv)]:
	for guide in find_guides(genome_seq):
		exclude_list.add(guide)

outfile = input("Please enter name of temporary buffer file for exlcude sequences): ")
out = open(outfile, "wt")

for item in exclude_list:
	out.write(f"{item}\n")

out.close()

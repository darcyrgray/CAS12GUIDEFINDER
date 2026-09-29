# CAS12GUIDEFINDER
Python program to find CAS12 guides from bacterial genomes for microbiology applications

Dependancy on Biopython repository.

Command line python program to help identify 25 base pair guide crRNA (including 4bp TTTV PAM) sequences within a target genome (can be entered as ASCII text file or FASTA).
Second function to exclude presence of these guides in other selected genomes.
This exclusion loop includes checking for specificity with 1 bp wobble - however run time with this is much higher so feel free to deactivate if needed.

All guides should still be checked with second database for specificty - I have used NCBI Blast to good effect.

Aimed at finding guides in microbial genomes.

Please note - author is a microbiology trainee trying to learn to implement computational biology tools including developing tools with BASH and Python. This script is not the most optimal search algorithm and may not be fit for all uses. 


Update 29/09/2026 - CAS12GF_rev.py is updated program - aimed to be combined with exclude_build.py --> uses more memory to build hash-table for faster look up for guide exclusion (using python sets instead of interating over whole genome/using lists).

check_SA.py and check_mec.py were uploaded for reference - used to check guide presence in other reference genomes.

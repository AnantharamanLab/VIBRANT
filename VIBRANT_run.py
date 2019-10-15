#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison, 2019

# VIBRANT v1.0.0
# Virus Identification By iteRative ANnoTation

# Usage: $ python3 VIBRANT_run.py -i <input_file> [options]

############################### Imports  #######################################
import warnings
warnings.filterwarnings("ignore")
import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
from datetime import date
import datetime
import argparse
import subprocess
import math
import time
import logging
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

vibrant_path = str(os.path.dirname(os.path.abspath(__file__)))

start = time.time()
round((time.time() - start)/60,1)

start_time = str(datetime.datetime.now().time()).rsplit(".",1)[0]

vibrant = argparse.ArgumentParser(description='Usage: python3 VIBRANT_run.py -i <input_file> [options]. VIBRANT identifies bacterial and archaeal viruses (phages) from assembled metagenomic scaffolds or whole genomes, including the excision of integrated proviruses. VIBRANT also performs curation of identified viral scaffolds, estimation of viral genome completeness and analysis of viral metabolic capabilities.')
vibrant.add_argument('--version', action='version', version='VIBRANT v1.0.0')

####### Required
vibrant.add_argument('-i', type=str, nargs=1, required=True, help='input fasta file')

####### Optional
vibrant.add_argument('-f', type=str, nargs=1, default='nucl', choices=['prot','nucl'], help='format of input [default="nucl"]')
vibrant.add_argument('-t', type=str, nargs=1, default='1', help='number of parallel VIBRANT runs, each occupies 1 CPU [default=1, max of 1 CPU per scaffold]')
vibrant.add_argument('-l', type=str, nargs=1, default='3000', help='length in basepairs to limit input sequences [default=3000, can increase but not decrease]')
vibrant.add_argument('-o', type=str, nargs=1, default='4', help='number of ORFs per scaffold to limit input sequences [default=4, can increase but not decrease]')
vibrant.add_argument('-virome', action='store_true', help='use this setting if dataset is known to be comprised mainly of viruses. More sensitive to viruses, less sensitive to false identifications [default=off]')
vibrant.add_argument('-no_plot', action='store_true', help='suppress the generation of summary plots [default=off]')
vibrant.add_argument('-k', type=str, nargs=1, default=str(vibrant_path) + '/databases/KEGG_profiles_prokaryotes.HMM', help='path to KEGG HMMs (if moved from default location)')
vibrant.add_argument('-p', type=str, nargs=1, default=str(vibrant_path) + '/databases/Pfam-A_v32.HMM', help='path to Pfam HMMs (if moved from default location)')
vibrant.add_argument('-v', type=str, nargs=1, default=str(vibrant_path) + '/databases/VOGDB94_phage.HMM', help='path to VOG HMMs (if moved from default location)')
vibrant.add_argument('-e', type=str, nargs=1, default=str(vibrant_path) + '/databases/Pfam-A_plasmid_v32.HMM', help='path to plasmid HMMs (if moved from default location)')
vibrant.add_argument('-a', type=str, nargs=1, default=str(vibrant_path) + '/databases/Pfam-A_phage_v32.HMM', help='path to viral-subset Pfam HMMs (if moved from default location)')
vibrant.add_argument('-c', type=str, nargs=1, default=str(vibrant_path) + '/files/VIBRANT_categories.tsv', help='path to VIBRANT categories file (if moved from default location)')
vibrant.add_argument('-n', type=str, nargs=1, default=str(vibrant_path) + '/files/VIBRANT_names.tsv', help='path to VIBRANT annotation to name file (if moved from default location)')
vibrant.add_argument('-s', type=str, nargs=1, default=str(vibrant_path) + '/files/VIBRANT_KEGG_pathways_summary.tsv', help='path to VIBRANT summary of KEGG metabolism file (if moved from default location)')
vibrant.add_argument('-m', type=str, nargs=1, default=str(vibrant_path) + '/files/VIBRANT_machine_model.sav', help='path to VIBRANT neural network machine learning model (if moved from default location)')
vibrant.add_argument('-g', type=str, nargs=1, default=str(vibrant_path) + '/files/VIBRANT_AMGs.tsv', help='path to VIBRANT AMGs file (if moved from default location)')
args = vibrant.parse_args()

if type(args.f) == str:
	format = args.f
else:
	format = "".join(args.f[0])
if args.virome == False:
    virome = ''
elif args.virome == True:
    virome = ' -virome'
suppress = args.no_plot
if type(args.k) == str:
	kegg_hmm = args.k
else:
	kegg_hmm = "".join(args.k[0])
if type(args.p) == str:
	pfam_hmm = args.p
else:
	pfam_hmm = "".join(args.p[0])
if type(args.v) == str:
	vog_hmm = args.v
else:
	vog_hmm = "".join(args.v[0])
if type(args.e) == str:
	plasmid_hmm = args.e
else:
	plasmid_hmm = "".join(args.e[0])
if type(args.a) == str:
	vpfam_hmm = args.a
else:
	vpfam_hmm = "".join(args.a[0])
if args.l == '3000':
	lim_low = args.l
if args.l != '3000':
	lim_low = args.l[0]
if args.t == '1':
	threads = int(args.t)
if args.t != '1':
	threads = int(args.t[0])
if args.o == '4':
	orf_low = args.o
if args.o != '4':
	orf_low = args.o[0]
if type(args.c) == str:
	categories = args.c
else:
	categories = "".join(args.c[0])
if type(args.n) == str:
	annotation_names = args.n
else:
	annotation_names = "".join(args.n[0])
if type(args.m) == str:
	model = args.m
else:
	model = "".join(args.m[0])
if type(args.g) == str:
	AMG_list = args.g
else:
	AMG_list = "".join(args.g[0])
if type(args.s) == str:
	summary_list = args.s
else:
	summary_list = "".join(args.s[0])
try:
    filename = str(args.i[0].rsplit('/')[-1])
    filepath = str(args.i[0].rsplit('/')[0])
except Exception:
    filename = str(args.i[0])

#######
if int(orf_low) < 4 or int(lim_low) < 3000:
	print("VIBRANT error: minimum ORFs is 4 and minimum sequence length is 3kb. These variables can only increase. Exiting." + "\n")
	exit()
command = ' -f ' + str(format) + ' -k ' + str(kegg_hmm) + ' -p ' + str(pfam_hmm) + ' -v ' + str(vog_hmm) + ' -e ' + str(plasmid_hmm) + ' -a ' + str(vpfam_hmm) + ' -c ' + str(categories) + ' -n ' + str(annotation_names) + ' -m ' + str(model) + ' -l ' + str(lim_low) + ' -o ' + str(orf_low) + str(virome) + ' -g ' + str(AMG_list)
base = str(filename.rsplit(".",1)[0])
logging.basicConfig(filename='VIBRANT_log_' + base + '.log', level=logging.INFO, format='%(message)s')
subprocess.run("echo '' > " + base + "_four-orf-count.txt", shell=True)

#######
if format == "prot":
    input_sequences = []
    with open(str(args.i[0]),'r') as infile:
        with open(str(base)+'.parallel-runs.txt', 'w') as outfile:
            for name, seq in SimpleFastaParser(infile):
                if " # " not in name:
                    subprocess.run('rm ' + str(base)+'.parallel-runs* 2>/dev/null', shell=True)
                    subprocess.run('rm ' + str(base)+'_four-orf-count.txt 2>/dev/null', shell=True)
                    subprocess.run('rm VIBRANT_log_' + str(base) + '.log 2>/dev/null', shell=True)
                    print("\n")
                    print("VIBRANT Error: either the wrong format (-f) was set or proteins are not in Prodigal format. Check input file (-i) and format (-f).")
                    print("See example data file 'Microviridae_ctch16_MH552510.2.faa' for Prodigal format.")
                    print("\n")
                    exit()
                protein = name.split(" # ",1)[0]
                genome = protein.rsplit("_",1)[0]
                input_sequences.append(genome)
                outfile.write(str(genome) + '\n')

    subprocess.run('cat ' + str(base)+'.parallel-runs.txt | uniq -d > ' + str(base)+'.parallel-runs.temp', shell=True)
    listing = 'wc -l ' + str(base)+'.parallel-runs.temp'
    listing_shell = subprocess.check_output(listing, shell=True)
    listing1 = str(listing_shell.strip()).split("'")[1]
    sequences = str(listing1.split("\\n")[0].split(" ",1)[0])

if format == "nucl":
    sequences = 0
    with open(str(args.i[0]),'r') as infile:
        with open(str(base)+'.parallel-runs.temp', 'w') as outfile:
            for name, seq in SimpleFastaParser(infile):
                outfile.write(str(name) + '\n')
                sequences += 1

shuffle = 'sort -R ' + str(base)+'.parallel-runs.temp > ' + str(base)+'.parallel-runs.txt'
s2 = subprocess.Popen(shuffle, shell=True)
s2.wait()
subprocess.Popen('rm ' + base + '.parallel-runs.temp', shell=True)
time.sleep(0.1)

lines = math.ceil(int(sequences)/int(threads))
parallel = 'split -l ' + str(lines) + ' ' + str(base)+'.parallel-runs.txt' + ' ' + str(base)+'.parallel-runs_'
subprocess.run(parallel, shell=True)
time.sleep(0.1)

move = 'for file in ' + str(base) + '.parallel-runs_*; do mv $file ${file}.txt; done'
s1 = subprocess.Popen(move, shell=True)
s1.wait()
time.sleep(0.1)

parallels = 'ls ' + str(base) + '.parallel-runs_*.txt'
parallels_shell = subprocess.check_output(parallels, shell=True)
parallels1 = str(parallels_shell.strip()).split("'")[1]
parallels2 = parallels1.split("\\n")

n = 0
p_list = []
while n < len(parallels2):
	if format == "nucl":
		execute = str(vibrant_path) + '/scripts/VIBRANT_extract_nucleotide.py ' + str(parallels2[n]) + " " + str(parallels2[n]).rsplit(".",1)[0] + ".fna " + str(args.i[0]) + " " + str(lim_low)
	elif format == "prot":
		execute = str(vibrant_path) + '/scripts/VIBRANT_extract_protein.py ' + str(parallels2[n]) + " " + str(parallels2[n]).rsplit(".",1)[0] + ".faa " + str(args.i[0])
	p = subprocess.Popen(execute, shell=True)
	p_list.append(p)
	n += 1

for item in p_list:
    item.wait()
time.sleep(0.1)

if format == "nucl":
	try:
		listing = 'ls ' + str(base) + '.parallel-runs_*.fna 2>/dev/null'
		listing_shell = subprocess.check_output(listing, shell=True)
		listing1 = str(listing_shell.strip()).split("'")[1]
		listing2 = listing1.split("\\n")
	except Exception:
		subprocess.run('rm ' + str(base)+'.parallel-runs* 2>/dev/null', shell=True)
		subprocess.run('rm ' + str(base)+'_four-orf-count.txt 2>/dev/null', shell=True)
		subprocess.run('rm VIBRANT_log_' + str(base) + '.log 2>/dev/null', shell=True)
		print("\n")
		print('No phages found. There were either no scaffolds at least ' + str(lim_low) + 'bp (set minimum size) or the input file (-i) format (-f) was not FASTA nucleotide. Exiting.')
		print("\n")
		exit()

if format == "prot":
	try:
		listing = 'ls ' + str(base) + '.parallel-runs_*.faa  2>/dev/null'
		listing_shell = subprocess.check_output(listing, shell=True)
		listing1 = str(listing_shell.strip()).split("'")[1]
		listing2 = listing1.split("\\n")
	except Exception:
		subprocess.run('rm ' + str(base)+'.parallel-runs* 2>/dev/null', shell=True)
		subprocess.run('rm ' + str(base)+'_four-orf-count.txt 2>/dev/null', shell=True)
		subprocess.run('rm VIBRANT_log_' + str(base) + '.log 2>/dev/null', shell=True)
		print("\n")
		print("No phages found. There were either no scaffolds at least ' + str(orf_low) + ' ORFs (set minimum size) or VIBRANT never found any proteins. Verify input file (-i) format (-f) is FASTA protein. Exiting.")
		exit()

remove = 'rm ' + str(base) + '.parallel-runs_*.txt'
subprocess.run(remove, shell=True)
time.sleep(0.1)

n = 0
p_list = []
while n < len(listing2):
    execute = str(vibrant_path) + '/scripts/VIBRANT_annotation.py -i ' + str(listing2[n]) + str(command)
    p = subprocess.Popen(execute, shell=True)
    p_list.append(p)
    n += 1

for item in p_list:
    item.wait()
time.sleep(0.1)

try:
	grep = "grep -c '1' " + base + "_four-orf-count.txt 2>/dev/null"
	grep_shell = subprocess.check_output(grep, shell=True)
	grep_out_count = str(grep_shell.strip()).split("'")[1]
	time.sleep(0.1)
except Exception:
	print("\n")
	print('No phages found. There were no scaffolds at least ' + str(orf_low) + ' ORFs (set minimum size). Exiting.')
	print("\n")
	subprocess.run('rm ' + str(base)+'.parallel-runs* 2>/dev/null', shell=True)
	subprocess.run('rm ' + str(base)+'_four-orf-count.txt 2>/dev/null', shell=True)
	subprocess.run('rm VIBRANT_log_' + str(base) + '.log 2>/dev/null', shell=True)
	exit()

if format == "nucl":
    cat_combo_ffn = 'cat ' + base + '*parallel-runs*.phages_combined.ffn 1> ' + str(base) + '.phages_combined.ffn 2>/dev/null'
    subprocess.Popen(cat_combo_ffn, shell=True)
    cat_combo_fna = 'cat ' + base + '*parallel-runs*.phages_combined.fna 1> ' + str(base) + '.phages_combined.fna 2>/dev/null'
    subprocess.Popen(cat_combo_fna, shell=True)
    cat_lytic_ffn = 'cat ' + base + '*parallel-runs*.phages_lytic.ffn 1> ' + str(base) + '.phages_lytic.ffn 2>/dev/null'
    subprocess.Popen(cat_lytic_ffn, shell=True)
    cat_lytic_fna = 'cat ' + base + '*parallel-runs*.phages_lytic.fna 1> ' + str(base) + '.phages_lytic.fna 2>/dev/null'
    r1 = subprocess.Popen(cat_lytic_fna, shell=True)
    cat_lysogenic_ffn = 'cat ' + base + '*parallel-runs*.phages_lysogenic.ffn 1> ' + str(base) + '.phages_lysogenic.ffn 2>/dev/null'
    subprocess.Popen(cat_lysogenic_ffn, shell=True)
    cat_lysogenic_fna = 'cat ' + base + '*parallel-runs*.phages_lysogenic.fna 1> ' + str(base) + '.phages_lysogenic.fna 2>/dev/null'
    r2 = subprocess.Popen(cat_lysogenic_fna, shell=True)
    r1.wait()
    r2.wait()

cat_combo_faa = 'cat ' + base + '*parallel-runs*.phages_combined.faa 1> ' + str(base) + '.phages_combined.faa 2>/dev/null'
r3 = subprocess.Popen(cat_combo_faa, shell=True)
cat_lytic_faa = 'cat ' + base + '*parallel-runs*.phages_lytic.faa 1> ' + str(base) + '.phages_lytic.faa 2>/dev/null'
r4 = subprocess.Popen(cat_lytic_faa, shell=True)
cat_lysogenic_faa = 'cat ' + base + '*parallel-runs*.phages_lysogenic.faa 1> ' + str(base) + '.phages_lysogenic.faa 2>/dev/null'
r5 = subprocess.Popen(cat_lysogenic_faa, shell=True)
r3.wait()
r4.wait()
r5.wait()
time.sleep(0.1)

if format == "nucl":
    with open(str(base) + '.phages_lysogenic.fna', 'r') as lysogenic_phages:
        with open(str(base) + '.phages_lytic.fna', 'r') as lytic_phages:
            with open('VIBRANT_complete_circular_' + base + '.tsv', 'w') as circular_phages:
                for name, seq in SimpleFastaParser(lysogenic_phages):
                    if "fragment" not in name:
                        seq = list(seq)
                        seq_start = "".join(seq[0:20])
                        n = int(len(seq))-1000
                        while n <= len(seq):
                            seq_end = "".join(seq[n:n+20])
                            if str(seq_start) == str(seq_end):
                                circular_phages.write(str(name) + '\t' + 'lysogenic' + '\t' + 'complete circular' + '\n')
                                n += 1000
                            else:
                                n += 1
                for item, sequence in SimpleFastaParser(lytic_phages):
                    sequence = list(sequence)
                    sequence_start = "".join(sequence[0:20])
                    n = int(len(sequence))-1000
                    while n <= len(sequence):
                        sequence_end = "".join(sequence[n:n+20])
                        if str(sequence_start) == str(sequence_end):
                            circular_phages.write(str(item) + '\t' + 'lytic' + '\t' + 'complete circular' + '\n')
                            n += 1000
                        else:
                            n += 1

f1 = subprocess.Popen('rm ' + base + '*parallel-runs*phages_*.faa 2>/dev/null', shell=True)
subprocess.Popen('rm ' + base + '*parallel-runs*phages_*.ffn 2>/dev/null', shell=True)
subprocess.Popen('rm ' + base + '*parallel-runs*phages_*.fna 2>/dev/null', shell=True)
f1.wait()
cat_faa = 'cat ' + base + '.parallel-runs_*.faa > ' + base + '.faa 2>/dev/null'
subprocess.Popen(cat_faa, shell=True)
cat_ffn = 'cat ' + base + '.parallel-runs_*.ffn > ' + base + '.ffn 2>/dev/null'
subprocess.Popen(cat_ffn, shell=True)

cat_annotations = 'cat VIBRANT_annotations*' + base + '* > temp_VIBRANT_annotations_' + str(base) + '.txt 2>/dev/null'
p6 = subprocess.Popen(cat_annotations, shell=True)
cat_results = 'cat VIBRANT_results*' + base + '* > temp_VIBRANT_results_' + str(base) + '.txt 2>/dev/null'
p7 = subprocess.Popen(cat_results, shell=True)
cat_amg = 'cat VIBRANT_AMGs*' + base + '* 1> temp_VIBRANT_AMGs_' + str(base) + '.txt 2>/dev/null'
p10 = subprocess.Popen(cat_amg, shell=True)
cat_gb = 'cat VIBRANT_genbank_table*' + base + '* 1> temp_VIBRANT_genbank_table_' + str(base) + '.txt 2>/dev/null'
p11 = subprocess.Popen(cat_gb, shell=True)
cat_machine = 'cat temp_VIBRANT_machine.*' + base + '* 1> temp2_VIBRANT_machine_' + str(base) + '.txt 2>/dev/null'
p8 = subprocess.Popen(cat_machine, shell=True)
rm_unmod = 'rm unmodified_VIBRANT_results*' + base + '* 2>/dev/null'
p9 = subprocess.Popen(rm_unmod, shell=True)
cat_quality = 'cat ' + base + '*genome_quality.out 1> temp_VIBRANT_genome_quality_' + base + '.tsv 2>/dev/null'
y11 = subprocess.Popen(cat_quality, shell=True)
cat_names = 'cat ' + base + '.parallel-runs_*.phages_combined.out 1> ' + base + '.phages_combined.txt 2>/dev/null'
y12 = subprocess.Popen(cat_names, shell=True)
y11.wait()
y12.wait()
p6.wait()
p7.wait()
p8.wait()
p9.wait()
p10.wait()
p11.wait()
time.sleep(0.1)
header_quality = 'echo "scaffold\ttype\tQuality" | cat - ' + 'temp_VIBRANT_genome_quality_' + base + '.tsv 1> ' + 'VIBRANT_genome_quality_' + base + '.tsv 2>/dev/null'
p7 = subprocess.Popen(header_quality, shell=True)
p7.wait()
time.sleep(0.1)

if format == "nucl":
	subprocess.run('mv VIBRANT_genome_quality_' + base + '.tsv temp_VIBRANT_genome_quality_' + base + '.tsv', shell=True)
	subprocess.run('cat temp_VIBRANT_genome_quality_' + base + '.tsv VIBRANT_complete_circular_' + base + '.tsv > VIBRANT_genome_quality_' + base + '.tsv', shell=True)

subprocess.Popen('rm VIBRANT_annotations*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm VIBRANT_results*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm VIBRANT_AMG*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm VIBRANT_genbank_table*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm temp_VIBRANT_machine.*' + base + '* 2>/dev/null', shell=True)

header_anno = 'echo "protein\tscaffold\tKO\tAMG\tKO name\tKO evalue\tKO score\tKO v-score\tPfam\tPfam name\tPfam evalue\tPfam score\tPfam v-score\tVOG\tVOG name\tVOG evalue\tVOG score\tVOG v-score" | cat - ' + 'temp_VIBRANT_annotations_' + str(base) + '.txt > ' + 'VIBRANT_annotations_' + str(base) + '.tsv'
p1 = subprocess.Popen(header_anno, shell=True)
header_results ='echo "scaffold\ttotal genes\tall KEGG\tKEGG v-score\tall Pfam\tPfam v-score\tall VOG\tVOG v-score\tKEGG int-rep\tKEGG zero\tPfam int-rep\tPfam zero\tVOG redoxin\tVOG rec-tran\tVOG int\tVOG RnR\tVOG DNA\tKEGG restriction check\tKEGG toxin check\tVOG special\tannotation check\tp_v check\tp_k check\tk_v check\tk check\tp check\tv check\th check" | cat - ' + 'temp_VIBRANT_results_' + str(base) + '.txt > ' + 'VIBRANT_summary_results_' + str(base) + '.tsv'
p2 = subprocess.Popen(header_results, shell=True)
header_amg = 'echo "protein\tscaffold\tAMG KO\tAMG KO name\tPfam\tPfam name" | cat - ' + 'temp_VIBRANT_AMGs_' + str(base) + '.txt 1> ' + 'VIBRANT_AMG_individuals_' + str(base) + '.tsv 2>/dev/null'
p4 = subprocess.Popen(header_amg, shell=True)
header_gb = 'echo "protein\tscaffold\taccession\tname" | cat - ' + 'temp_VIBRANT_genbank_table_' + str(base) + '.txt 1> ' + 'VIBRANT_genbank_table_' + str(base) + '.tsv 2>/dev/null'
p6 = subprocess.Popen(header_gb, shell=True)
for1 = 'for f in ' + base + '*parse.txt; do sed "1d" $f > ${f%.*}.out; done'
p3 = subprocess.Popen(for1, shell=True)
header_machine = 'echo "scaffold\tprediction" | cat - ' + 'temp2_VIBRANT_machine_' + str(base) + '.txt 1> ' + 'VIBRANT_machine_' + str(base) + '.tsv 2>/dev/null'
p5 = subprocess.Popen(header_machine, shell=True)
if format == "nucl":
    for2 = 'for f in ' + base + '*parallel-runs*temp; do sed "1d" $f 1> ${f%.*}.temp.gff 2>/dev/null; done'
    subprocess.run(for2, shell=True)
p1.wait()
p2.wait()
p3.wait()
p4.wait()
p5.wait()
p6.wait()
time.sleep(0.1)

subprocess.Popen('rm temp_VIBRANT_annotations*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm temp_VIBRANT_results*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm temp2_VIBRANT_machine*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm temp_VIBRANT_AMGs_*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm temp_VIBRANT_genbank_table_*' + base + '* 2>/dev/null', shell=True)
subprocess.Popen('rm temp_VIBRANT_genome_quality_' + base + '.tsv 2>/dev/null', shell=True)
subprocess.Popen(['mkdir', 'VIBRANT_HMM_tables_unformatted_'+base])
subprocess.Popen(['mkdir', 'VIBRANT_HMM_tables_parsed_'+base])
subprocess.Popen(['mkdir', 'VIBRANT_phages_'+base])
subprocess.Popen(['mkdir', 'VIBRANT_results_'+base])
kegg_cat = 'cat ' + base + '*KEGG.hmmtbl.parse.out 1> ' + base + '.KEGG_hmmtbl_parse.txt 2>/dev/null'
x1 = subprocess.Popen(kegg_cat, shell=True)
kegg_table = 'cat ' + base + '*KEGG.hmmtbl 1> ' + base + '_unformatted_KEGG.hmmtbl 2>/dev/null'
y1 = subprocess.Popen(kegg_table, shell=True)
pfam_cat = 'cat ' + base + '*Pfam.hmmtbl.parse.out 1> ' + base + '.Pfam_hmmtbl_parse.txt 2>/dev/null'
x2 = subprocess.Popen(pfam_cat, shell=True)
pfam_table = 'cat ' + base + '*Pfam.hmmtbl 1> ' + base + '_unformatted_Pfam.hmmtbl 2>/dev/null'
y2 = subprocess.Popen(pfam_table, shell=True)
vog_cat = 'cat ' + base + '*VOG.hmmtbl.parse.out 1> ' + base + '.VOG_hmmtbl_parse.txt 2>/dev/null'
x3 = subprocess.Popen(vog_cat, shell=True)
vog_table = 'cat ' + base + '*VOG.hmmtbl 1> ' + base + '_unformatted_VOG.hmmtbl 2>/dev/null'
y3 = subprocess.Popen(vog_table, shell=True)
vpfam_cat = 'cat ' + base + '*vpfam.hmmtbl.parse.out 1> ' + base + '.vpfam_hmmtbl_parse.txt 2>/dev/null'
vpfam_table = 'cat ' + base + '*vpfam.hmmtbl 1> ' + base + '_unformatted_vpfam.hmmtbl 2>/dev/null'
x4 = subprocess.Popen(vpfam_cat, shell=True)
y4 = subprocess.Popen(vpfam_table, shell=True)
plasmid_cat = 'cat ' + base + '*plasmid.hmmtbl.parse.out 1> ' + base + '.plasmid_hmmtbl_parse.txt 2>/dev/null'
plasmid_table = 'cat ' + base + '*plasmid.hmmtbl 1> ' + base + '_unformatted_plasmid.hmmtbl 2>/dev/null'
x5 = subprocess.Popen(plasmid_cat, shell=True)
y5 = subprocess.Popen(plasmid_table, shell=True)
x1.wait()
x2.wait()
x3.wait()
x4.wait()
y1.wait()
y2.wait()
y3.wait()
y4.wait()
y5.wait()
x5.wait()
time.sleep(0.1)

move_parse = 'mv ' + base + '*_hmmtbl_parse.txt ' + 'VIBRANT_HMM_tables_parsed_'+base + ' 2>/dev/null'
x6 = subprocess.Popen(move_parse, shell=True)
move_table = 'mv ' + base + '*_unformatted_*.hmmtbl ' + 'VIBRANT_HMM_tables_unformatted_'+base + ' 2>/dev/null'
y6 = subprocess.Popen(move_table, shell=True)
x6.wait()
y6.wait()
time.sleep(0.1)

for2 = 'for f in VIBRANT_HMM_tables_parsed_' + base + '/*.txt; do echo "protein\tid\tevalue\tscore" | cat - $f > ${f%.*}.tsv; done'
z1 = subprocess.Popen(for2, shell=True)
z1.wait()
rm_parsed = 'rm VIBRANT_HMM_tables_parsed_' + base + '/*.txt 2>/dev/null'
z2 = subprocess.Popen(rm_parsed, shell=True)
z2.wait()
time.sleep(0.1)

rm_faa = 'rm ' + base + '*parallel-runs*faa 2>/dev/null'
rm_fna = 'rm ' + base + '*parallel-runs*fna 2>/dev/null'
rm_ffn = 'rm ' + base + '*parallel-runs*ffn 2>/dev/null'
rm_list = 'rm ' + base + '*parallel-runs_*txt 2>/dev/null'
rm_names = 'rm ' + base + '*parallel-runs_*out 2>/dev/null'
rm_temp = 'rm ' + base + '*parallel-runs_*temp 2>/dev/null'
subprocess.Popen(rm_faa, shell=True)
subprocess.Popen(rm_fna, shell=True)
subprocess.Popen(rm_ffn, shell=True)
subprocess.Popen(rm_list, shell=True)
subprocess.Popen(rm_temp, shell=True)
subprocess.Popen(rm_names, shell=True)

cat1_faa = "cat " + base + ".faa | sed 's/$~&/ /g' 1> " + base + ".temp.faa 2>/dev/null"
cat2_faa = "cat " + base + ".temp.faa | sed 's/^@%/\"/g' 1> " + base + ".faa 2>/dev/null"
rm1 = 'rm ' + base + '.faa 2>/dev/null'
subprocess.run(cat1_faa, shell=True)
subprocess.run(rm1, shell=True)
subprocess.run(cat2_faa, shell=True)
subprocess.run("rm " + base + ".temp.faa 2>/dev/null", shell=True)
time.sleep(0.1)

if format == "nucl":
    cat1_ffn = "cat " + base + ".ffn | sed 's/$~&/ /g' 1> " + base + ".temp.ffn 2>/dev/null"
    cat2_ffn = "cat " + base + ".temp.ffn | sed 's/^@%/\"/g' 1> " + base + ".ffn 2>/dev/null"
    cat_gff = "cat " + base + "*temp.gff 1> " + base + ".all.gff 2>/dev/null"
    c2 = subprocess.Popen(cat1_ffn, shell=True)
    c7 = subprocess.Popen(cat_gff, shell=True)
    c2.wait()
    c7.wait()
    time.sleep(0.1)
    cat1_gff = "cat " + base + ".all.gff | sed 's/$~&/ /g' 1> " + base + ".temp.gff 2>/dev/null"
    cat2_gff = "cat " + base + ".temp.gff | sed 's/^@%/\"/g' 1> " + base + ".all.gff 2>/dev/null"
    c10 = subprocess.Popen(cat1_gff, shell=True)
    c11 = subprocess.Popen(cat2_gff, shell=True)
    c10.wait()
    c11.wait()
    time.sleep(0.1)
    rm2 = 'rm ' + base + '.ffn 2>/dev/null'
    rm3 = 'rm ' + base + '*temp.gff'
    c6 = subprocess.Popen(rm2, shell=True)
    c8 = subprocess.Popen(rm3, shell=True)
    c6.wait()
    c8.wait()
    c4 = subprocess.Popen(cat2_ffn, shell=True)
    header_gff = 'echo "##gff-version  3" | cat - ' + base + '.all.gff 1> ' + base + '.gff 2>/dev/null'
    c9 = subprocess.Popen(header_gff, shell=True)
    c4.wait()
    c9.wait()
    subprocess.Popen("rm " + base + ".temp.ffn 2>/dev/null", shell=True)
    subprocess.Popen("rm " + base + ".all.gff 2>/dev/null", shell=True)

with open('VIBRANT_summary_results_'+base+'.tsv') as results:
    out_phages = -1
    for line in results:
        out_phages += 1

remove_run1 = str('rm ' + base + '.parallel-runs*hmmtbl 2>/dev/null')
subprocess.Popen(remove_run1, shell=True)
remove_run2 = str('rm ' + base + '.parallel-runs*hmmtbl.parse* 2>/dev/null')
subprocess.Popen(remove_run2, shell=True)

genbank_database = []
with open(base + '.phages_combined.faa', 'r') as fasta:
	for name, seq in SimpleFastaParser(fasta):
		sequence = ''
		if len(seq) > 44:
			sequence += str(seq[0:44])+"\n"
			for group in range(44,len(seq),58):
				sequence += '                     ' + str(seq[group:group+58]) + '\n'
			sequence = sequence.rstrip("\n")
		elif len(seq) <= 44:
			sequence = seq
		protein = str(name.split("\t")[0])
		genome = str(protein.rsplit("_",1)[0])
		location = str(name.split("\t")[1])
		strand = str(name.split("\t")[2])
		accession = str(name.split("\t")[3])
		name = str(name.split("\t")[4])
		genbank_database.extend((protein,genome,location,strand,accession,name,str(sequence)))
	genbank_database.extend(('',''))

if format == "nucl":
	length_database = {}
	sequence_database = {}
	length_list = []
	with open(base + '.phages_combined.fna', 'r') as fasta:
		for name, seq in SimpleFastaParser(fasta):
			length_database.update({name:len(seq)})
			length_list.append(int(len(seq)))
			n = 0
			count = 0
			i = 61
			spaces_database = '               '
			sequence = '         1 '
			for group in range(0,len(seq),10):
				sequence += str(seq[group:group+10]) + ' '
				count += 1
				if count == 6:
					length = list(str(i))
					size = 9-len(length)
					spaces = str(spaces_database[0:size+1])
					sequence += '\n' + str(spaces) + str(i) + ' '
					i += 60
					count = 0
			sequence_database.update({name:sequence})

if format == "nucl":
	with open('VIBRANT_complete_circular_' + base + '.tsv', 'r') as circular:
		circular_database = circular.read().replace('\n','\t').split('\t')

	with open(base + '.phages_combined.temp.gbk', 'w') as genbank:
		n = 1
		check = False
		while n < len(genbank_database)-1:
			if check == False:
				start_site = int(str(genbank_database[n+1][1:-1]).split(".",1)[0])-1
				if str(genbank_database[n]) in circular_database:
					form = "circular"
				else:
					form = "linear"
				genbank.write('//\nLOCUS       ' + str(genbank_database[n]) + '                 ' + str(length_database[genbank_database[n]]) + ' bp    DNA     ' + str(form) + '   VRL ' + str(date.today()) + '\nDEFINITION  ' + str(genbank_database[n]) + '.\nFEATURES             Location/Qualifiers\n     source          /organism="' + str(genbank_database[n]) + '"\n')
				check = True
			if "fragment" in str(genbank_database[n]):
				if check == True:
					location_A = int(str(genbank_database[n+1][1:-1]).split("..")[0])-start_site
					location_B = int(str(genbank_database[n+1][1:-1]).split("..")[1])-start_site
					location = str(location_A) + ".." + str(location_B)
					if str(genbank_database[n+2]) == "-1":
						if str(genbank_database[n]) == str(genbank_database[n+7]):
							genbank.write('     CDS             complement(' + str(location) + ')' + '\n                     /locus_tag="' + str(genbank_database[n-1]) + '"\n                     /protein_id="' + str(genbank_database[n+3]) + '"\n                     /product="' + str(genbank_database[n+4]) + '"\n                     /translation="' + str(genbank_database[n+5]) + '"\n')
						if str(genbank_database[n]) != str(genbank_database[n+7]):
							genbank.write('     CDS             complement(' + str(location) + ')' + '\n                     /locus_tag="' + str(genbank_database[n-1]) + '"\n                     /protein_id="' + str(genbank_database[n+3]) + '"\n                     /product="' + str(genbank_database[n+4]) + '"\n                     /translation="' + str(genbank_database[n+5]) + '"\n')
							genbank.write('ORIGIN\n' + str(sequence_database[genbank_database[n]]) + '\n')
							check = False
							start_site = 0
					elif str(genbank_database[n+2]) == "1":
						if str(genbank_database[n]) == str(genbank_database[n+7]):
							genbank.write('     CDS             ' + str(location) + '\n                     /locus_tag="' + str(genbank_database[n-1]) + '"\n                     /protein_id="' + str(genbank_database[n+3]) + '"\n                     /product="' + str(genbank_database[n+4]) + '"\n                     /translation="' + str(genbank_database[n+5]) + '"\n')
						if str(genbank_database[n]) != str(genbank_database[n+7]):
							genbank.write('     CDS             ' + str(location) + '\n                     /locus_tag="' + str(genbank_database[n-1]) + '"\n                     /protein_id="' + str(genbank_database[n+3]) + '"\n                     /product="' + str(genbank_database[n+4]) + '"\n                     /translation="' + str(genbank_database[n+5]) + '"\n')
							genbank.write('ORIGIN\n' + str(sequence_database[genbank_database[n]]) + '\n')
							check = False
							start_site = 0
			elif "fragment" not in str(genbank_database[n]):
				if check == True:
					if str(genbank_database[n+2]) == "-1":
						if str(genbank_database[n]) == str(genbank_database[n+7]):
							genbank.write('     CDS             complement' + str(genbank_database[n+1]) + '\n                     /locus_tag="' + str(genbank_database[n-1]) + '"\n                     /protein_id="' + str(genbank_database[n+3]) + '"\n                     /product="' + str(genbank_database[n+4]) + '"\n                     /translation="' + str(genbank_database[n+5]) + '"\n')
						if str(genbank_database[n]) != str(genbank_database[n+7]):
							genbank.write('     CDS             complement' + str(genbank_database[n+1]) + '\n                     /locus_tag="' + str(genbank_database[n-1]) + '"\n                     /protein_id="' + str(genbank_database[n+3]) + '"\n                     /product="' + str(genbank_database[n+4]) + '"\n                     /translation="' + str(genbank_database[n+5]) + '"\n')
							genbank.write('ORIGIN\n' + str(sequence_database[genbank_database[n]]) + '\n')
							check = False
					elif str(genbank_database[n+2]) == "1":
						if str(genbank_database[n]) == str(genbank_database[n+7]):
							genbank.write('     CDS             ' + str(genbank_database[n+1][1:-1]) + '\n                     /locus_tag="' + str(genbank_database[n-1]) + '"\n                     /protein_id="' + str(genbank_database[n+3]) + '"\n                     /product="' + str(genbank_database[n+4]) + '"\n                     /translation="' + str(genbank_database[n+5]) + '"\n')
						if str(genbank_database[n]) != str(genbank_database[n+7]):
							genbank.write('     CDS             ' + str(genbank_database[n+1][1:-1]) + '\n                     /locus_tag="' + str(genbank_database[n-1]) + '"\n                     /protein_id="' + str(genbank_database[n+3]) + '"\n                     /product="' + str(genbank_database[n+4]) + '"\n                     /translation="' + str(genbank_database[n+5]) + '"\n')
							genbank.write('ORIGIN\n' + str(sequence_database[genbank_database[n]]) + '\n')
							check = False
			n += 7


with open(annotation_names, 'r') as names:
	names = names.read().replace('\n','\t').split('\t')
	names_dict = {}
	n = 0
	while n < len(names):
		names_dict.update({names[n]:names[n+1]})
		n += 2

with open(summary_list, 'r') as summary:
	with open('VIBRANT_AMG_individuals_' + str(base) + '.tsv', 'r') as annotations:
		with open('VIBRANT_AMG_counts_' + str(base) + '.tsv', 'w') as count_file:
			with open('VIBRANT_AMG_pathways_' + str(base) + '.tsv', 'w') as path_file:
				count_file.write('AMG count\tAMG KO\tAMG KO name\n')
				summary = summary.read().replace('\n','\t').split('\t')
				annotations = pd.read_csv(annotations, sep='\t', header=0)
				annotations = pd.DataFrame(annotations)
				temp_anno = list(annotations['AMG KO'])
				kegg_anno = []
				count_dict = {}
				for item in temp_anno:
					if str(item).startswith('K'):
						kegg_anno.append(item)
				kegg_set = list(set(kegg_anno))
				for item in kegg_set:
					count = kegg_anno.count(item)
					count_dict.update({str(item):str(count)})
					count_file.write(str(count) + "\t" + str(item) + "\t" + str(names_dict[item]) + "\n")

				n = 4
				sum_count = 0
				present = []
				path_file.write("KEGG Entry\tMetabolism\tPathway\tTotal AMGs\tPresent AMG KOs\n")
				while n < len(summary):
					KOs = summary[n+3].split("~")
					for item in KOs:
						if item in count_dict.keys():
							sum_count += int(count_dict[item])
							present.append(item)
					if sum_count > 0:
						path_file.write(str(summary[n]) + "\t" + str(summary[n+1]) + "\t" + str(summary[n+2]) + "\t" + str(sum_count) + "\t" + str(",".join(present)) + "\n")
					sum_count = 0
					present = []
					n += 4

sed_gb = 'sed "1d" ' + base + '.phages_combined.temp.gbk 1> ' + base + '.phages_combined.gbk 2>/dev/null'
y12 = subprocess.Popen(sed_gb, shell=True)
y12.wait()
subprocess.run('rm ' + base + '.phages_combined.temp.gbk 2>/dev/null',shell=True)



####### Visualization
if suppress == False and int(out_phages) > 0:
	plt.figure()
	carbohydrate = ['map00010', 'map00020', 'map00030', 'map00040', 'map00051', 'map00052', 'map00053', 'map00500', 'map00520', 'map00620', 'map00630', 'map00640', 'map00650', 'map00660', 'map00562']
	energy = ['map00190', 'map00195', 'map00196', 'map00710', 'map00720', 'map00680', 'map00910', 'map00920']
	lipid = ['map00061', 'map00062', 'map00071', 'map00072', 'map00073', 'map00100', 'map00120', 'map00121', 'map00140', 'map00561', 'map00564', 'map00565', 'map00600', 'map00590', 'map00591', 'map00592', 'map01040']
	nucleotide = ['map00230', 'map00240']
	amino = ['map00250', 'map00260', 'map00270', 'map00280', 'map00290', 'map00300', 'map00310', 'map00220', 'map00330', 'map00340', 'map00350', 'map00360', 'map00380', 'map00400', 'map00410', 'map00430', 'map00440', 'map00450', 'map00460', 'map00471', 'map00472', 'map00473', 'map00480']
	glycan = ['map00510', 'map00513', 'map00512', 'map00515', 'map00514', 'map00532', 'map00534', 'map00533', 'map00531', 'map00563', 'map00601', 'map00603', 'map00604', 'map00540', 'map00550', 'map00511', 'map00571', 'map00572']
	cofactor = ['map00730', 'map00740', 'map00750', 'map00760', 'map00770', 'map00780', 'map00785', 'map00790', 'map00670', 'map00830', 'map00860', 'map00130']
	terpenoid = ['map00900', 'map00902', 'map00909', 'map00904', 'map00906', 'map00905', 'map00981', 'map00908', 'map00903', 'map00281', 'map01052', 'map00522', 'map01051', 'map01059', 'map01056', 'map01057', 'map00253', 'map00523', 'map01054', 'map01053', 'map01055']
	secondary = ['map00940', 'map00945', 'map00941', 'map00944', 'map00942', 'map00943', 'map00901', 'map00403', 'map00950', 'map00960', 'map01058', 'map00232', 'map00965', 'map00966', 'map00402', 'map00311', 'map00332', 'map00261', 'map00331', 'map00521', 'map00524', 'map00525', 'map00231', 'map00401', 'map00404', 'map00405', 'map00333', 'map00254', 'map00998', 'map00999']
	aromatic = ['map00362', 'map00627', 'map00364', 'map00625', 'map00361', 'map00623', 'map00622', 'map00633', 'map00642', 'map00643', 'map00791', 'map00930', 'map00363', 'map00621', 'map00626', 'map00624', 'map00365', 'map00984', 'map00980', 'map00982', 'map00983']
	relay = ['map04122']

	with open('VIBRANT_AMG_pathways_' + str(base) + '.tsv', 'r') as pathways:
		pathways = pathways.read().replace('\n','\t').split('\t')
		if len(pathways) > 6:
			carbohydrate_count = 0
			energy_count = 0
			lipid_count = 0
			nucleotide_count = 0
			amino_count = 0
			glycan_count = 0
			cofactor_count = 0
			terpenoid_count = 0
			secondary_count = 0
			aromatic_count = 0
			relay_count = 0

			n = 5
			while n < len(pathways):
				if str(pathways[n]) in carbohydrate:
					carbohydrate_count += int(pathways[n+3])
				elif str(pathways[n]) in energy:
					energy_count += int(pathways[n+3])
				elif str(pathways[n]) in lipid:
					lipid_count += int(pathways[n+3])
				elif str(pathways[n]) in nucleotide:
					nucleotide_count += int(pathways[n+3])
				elif str(pathways[n]) in amino:
					amino_count += int(pathways[n+3])
				elif str(pathways[n]) in glycan:
					glycan_count += int(pathways[n+3])
				elif str(pathways[n]) in cofactor:
					cofactor_count += int(pathways[n+3])
				elif str(pathways[n]) in terpenoid:
					terpenoid_count += int(pathways[n+3])
				elif str(pathways[n]) in secondary:
					secondary_count += int(pathways[n+3])
				elif str(pathways[n]) in aromatic:
					aromatic_count += int(pathways[n+3])
				elif str(pathways[n]) in relay:
					relay_count += int(pathways[n+3])
				n += 5

			dataframe = pd.DataFrame({"AMG Metabolism Category": ['carbohydrates', 'energy metabolism', 'lipids', 'nucleotides', 'amino acids', 'glycan', 'cofactors/vitamins', 'terpenoids/polyketides', 'secondary metabolites', 'aromatic compounds', 'sulfur relay'], "Number of AMGs": [carbohydrate_count, energy_count, lipid_count, nucleotide_count, amino_count, glycan_count, cofactor_count, terpenoid_count, secondary_count, aromatic_count, relay_count]})
			figure = sns.barplot(x="AMG Metabolism Category", y="Number of AMGs", data=dataframe, palette=sns.color_palette("pastel"), edgecolor=sns.color_palette("bright"))
			sns.set(style="white")
			plt.ylim(0, None)
			plt.xticks(rotation=90)
			plt.tight_layout()
			plt.savefig('VIBRANT_figure_pathways_' + str(base) + '.pdf', format='pdf')

if suppress == False and int(out_phages) > 0:
	plt.figure()
	with open('VIBRANT_genome_quality_' + base + '.tsv', 'r') as qual_file:
		data_frame = pd.read_csv(qual_file, sep='\t', header=0)
		data_frame = pd.DataFrame(data_frame)
		sns.set(style="white")
		qual_fig = sns.countplot(x="Quality", data=data_frame, palette=sns.color_palette("pastel"), edgecolor=sns.color_palette("bright"))
		plt.xticks(rotation=90)
		plt.tight_layout()
		plt.savefig('VIBRANT_figure_quality_' + str(base) + '.pdf', format='pdf')

if suppress == False and format == "nucl" and int(out_phages) > 0:
	plt.figure()
	size = int(len(length_list)*0.1+1)
	if len(length_list) <= 10:
		size = int(len(length_list))
	if len(length_list) > 10 and len(length_list) <= 50:
		size = int(len(length_list)/2)
	if size > 50:
		size = 50
	if size < 1:
		size = 1
	sns.set(style="white")
	figure = sns.distplot(length_list, bins=size, kde=False, rug=False, color='green')
	plt.xlabel("Scaffold size")
	plt.ylabel("Count")
	plt.tight_layout()
	plt.savefig('VIBRANT_figure_sizes_' + str(base) + '.pdf', format='pdf')

if suppress == False:
	if int(sequences) >= 10:
		plt.figure()
		proportions = [int(sequences),int(grep_out_count),int(out_phages)]
		names = ["total sequences\n"+str(sequences), "correct size\n"+str(grep_out_count), "total phages\n"+str(out_phages)]
		sns.set(style="white", font_scale=0.75)
		chart = sns.color_palette('Blues', 3)
		fig, figure = plt.subplots()
		for n, p in enumerate(proportions):
		    circle = plt.Circle((1, p), radius=p, facecolor=chart[n], edgecolor="black")
		    figure.add_artist(circle)
		    if n == 0:
		        figure.text(p+(p/3), p*1.25, names[n], color=chart[n], horizontalalignment='center')
		        save = float(p)
		    if n == 1:
		        figure.text(save+(save/3), save, names[n], color=chart[n], horizontalalignment='center')
		    if n == 2:
		        figure.text(save+(save/3), save*0.75, names[n], color=chart[n], horizontalalignment='center')
		figure.set_xlim(-int(sequences)-(save*0.1), int(sequences)+(save*0.1))
		figure.set_ylim(-(save*0.1), 2*int(sequences)+(save*0.1))
		figure.set_aspect('equal')
		plt.tight_layout()
		figure.set_xticks([])
		figure.set_yticks([])
		sns.despine(fig=None, top=True, right=True, left=True, bottom=True, offset=None, trim=True)
		plt.savefig('VIBRANT_figure_phages_' + str(base) + '.pdf', format='pdf')

	subprocess.run('mkdir VIBRANT_figures_'+base, shell=True)

subprocess.run("rm " + base + "_four-orf-count.txt 2>/dev/null", shell=True)
subprocess.run("rm " + base + "*parallel-runs* 2>/dev/null", shell=True)
move_phages = 'mv ' + base + '.phages_* ' + 'VIBRANT_phages_'+base + ' 2>/dev/null'
move_run4 = str('mv VIBRANT*' + base + '.tsv ' + 'VIBRANT_results_' +base + ' 2>/dev/null')
move_figs = str('mv VIBRANT_figure_*' + base + '.pdf ' + 'VIBRANT_figures_' +base + ' 2>/dev/null')
y8 = subprocess.Popen(move_run4, shell=True)
y7 = subprocess.Popen(move_phages, shell=True)
y9 = subprocess.Popen(move_figs, shell=True)
y7.wait()
y8.wait()
y9.wait()
time.sleep(0.1)
#
#
#
log_command = sys.argv
logging.info("Command:  " + str(" ".join(log_command)))
logging.info("Date:     " + str(date.today()))
logging.info("Start:    " + str(start_time))
logging.info("End:      " + str(datetime.datetime.now().time()).rsplit(".",1)[0])
logging.info("Runtime:  " + str(round((time.time() - float(start))/60,1)) + " minutes")
logging.info("Program:  VIBRANT v1.0.0")
logging.info("\n")
logging.info(str(sequences) + " scaffolds were read in.")
if format == "nucl":
	logging.info(str(grep_out_count) + " scaffolds met minimum requirements: at least " + str(lim_low) + "bp and " + str(orf_low) + " ORFs.")
if format == "prot":
	logging.info(str(grep_out_count) + " scaffolds met the minimum requirement: at least " + str(orf_low) + " ORFs.")
logging.info(str(out_phages) + " putative phages were identified.")
logging.info("\n")
logging.info('                                                               ##')
logging.info('                                                             ##  ##')
logging.info('                                                           ##      ##')
logging.info('######   ##  ##     ##     #######   ######    #####       ##      ##')
logging.info('##  ##   ##  ##   ##  ##   ##        ##       ##             ##  ##')
logging.info('######   ######   ######   ##  ###   ######    ###             ##')
logging.info('##       ##  ##   ##  ##   ##   ##   ##           ##           ##')
logging.info('##       ##  ##   ##  ##   #######   ######   #####            ##')
logging.info('                                                            #  ##  #')
logging.info('                                                           # # ## # #')
logging.info('                                                          #   #  #   #')
logging.info('                                                         #            #')
logging.info("\n")
#
#
#
subprocess.run('mkdir VIBRANT_' + base, shell=True)
if format == "nucl":
	subprocess.run('mv ' + base + '.faa VIBRANT_' + base +' 2>/dev/null', shell=True)
	subprocess.run('mv ' + base + '.ffn VIBRANT_' + base +' 2>/dev/null', shell=True)
	subprocess.run('mv ' + base + '.gff VIBRANT_' + base +' 2>/dev/null', shell=True)
subprocess.run('mv VIBRANT*' + base + ' VIBRANT_' + base +' 2>/dev/null', shell=True)
subprocess.run('mv VIBRANT*' + base + '.log VIBRANT_' + base +' 2>/dev/null', shell=True)
time.sleep(0.1)
#
#
#

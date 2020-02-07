#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison

# VIBRANT v1.1.0
# Virus Identification By iteRative ANnoTation
# Release date: Feb 7 2020

# Usage: see VIBRANT_run.py

############################### Imports  #######################################
import sys
import pandas as pd
import argparse
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import OrderedDict
import math
import numpy as np
from sklearn import preprocessing
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report,confusion_matrix
import pickle

############################### Set Arguments  #################################
vibrant = argparse.ArgumentParser(description='See main wrapper script: VIBRANT_run.py. This script performs the bulk of the work but is not callable on its own.')
vibrant.add_argument('--version', action='version', version='VIBRANT v1.1.0')

####### Required
vibrant.add_argument('-i', type=str, nargs=1, required=True, help='input fasta file')

####### Optional
vibrant.add_argument('-f', type=str, nargs=1, default='nucl', choices=['prot','nucl'], help='format of input [default="nucl"]')
vibrant.add_argument('-l', type=str, nargs=1, default='1000', help='length in basepairs to limit input sequences [default=1000, can increase but not decrease]')
vibrant.add_argument('-o', type=str, nargs=1, default='4', help='number of ORFs per scaffold to limit input sequences [default=4, can increase but not decrease]')
vibrant.add_argument('-virome', action='store_true', help='use this setting if dataset is known to be comprised mainly of viruses. More sensitive to viruses, less sensitive to false identifications [default=off]')
vibrant.add_argument('-d', type=str, nargs=1, help='path to "databases" directory that contains .HMM files (if moved from default location)')
vibrant.add_argument('-m', type=str, nargs=1, help='path to "files" directory that contains .tsv and model files (if moved from default location)')

####### Create variables
args = vibrant.parse_args()
input = str(args.i[0])
kegg_hmm = str(args.d[0]) + 'KEGG_profiles_prokaryotes.HMM'
pfam_hmm = str(args.d[0]) + 'Pfam-A_v32.HMM'
vog_hmm = str(args.d[0]) + 'VOGDB94_phage.HMM'
plasmid_hmm = str(args.d[0]) + 'Pfam-A_plasmid_v32.HMM'
vpfam_hmm = str(args.d[0]) + 'Pfam-A_phage_v32.HMM'
virome = args.virome
lim_low = int(args.l[0])
orf_low = int(args.o[0])
categories = str(args.m[0]) + 'VIBRANT_categories.tsv'
AMG_list = str(args.m[0]) + 'VIBRANT_AMGs.tsv'
annotation_names = str(args.m[0]) + 'VIBRANT_names.tsv'
model = str(args.m[0]) + 'VIBRANT_machine_model.sav'
format = args.f[0]

cpu = "1"
grep_shell = subprocess.check_output("hmmsearch -h", shell=True)
grep_out = str(grep_shell.strip()).split("\\n")
for item in grep_out:
	if "--cpu" in item and "[" in item:
		cpu = "0"

############################### Set input  #####################################
infile = str((input.rsplit('.',1)[:-1])[0])
base = str(infile.rsplit("/",1)[-1])
in_base = str(infile.rsplit(".",1)[0]).rsplit("/",1)[-1]
path = str(infile.rsplit("/",1)[:-1][0])+"/"

############################### Limit lengths ##################################
if format == "nucl":
	subprocess.call(['prodigal', '-i', input, '-a', infile+'.faa', '-d', infile+'.ffn', '-p', 'meta', '-f', 'gff', '-q', '-o', infile+'.temp'])

#######
if format == "prot":
	in_proteins = input
elif format == "nucl":
	in_proteins = infile+'.faa'

####### create list of non-redundant genomes
with open(str(in_proteins), 'r') as write_fasta:
	names = []
	basenames = []
	correct = []
	base_count = 0
	n = 0
	check_first = 0
	for name, seq in SimpleFastaParser(write_fasta):
		if format == "prot":
			name = (str(name).split(" # ")[0]+"_"+str(name).split(" # ")[3]).replace(" ", "$~&").replace('"','^@%')
		elif format == "nucl":
			name = (str(name).split(" # ")[0]+"_"+str(name).split(" # ")[3])
		basename = str(name.rsplit("_",2)[0])
		names.append(name)
		basenames.append(basename)
		if len(basenames) == 0:
			exit()
		if check_first != 0:
			if n <= len(basenames):
				if basenames[n] == basenames[n+1]:
					base_count += 1
				else:
					base_count += 1
					if base_count >= int(orf_low):
						correct.append(basenames[n])
					base_count = 0
		n += 1
		if check_first == 0:
			n -= 1
		check_first += 1

# Last contig is not added in previous section
if base_count >= int(orf_low)-1:
	correct.append(basenames[-1])

if len(correct) == 0:
	exit()

for item in correct:
	if item != "":
		subprocess.run("echo '1' >> " + str(path)+str(in_base) + "_four-orf-count.txt", shell=True)

####### write out genomes with correct number of ORFs
with open(str(in_proteins), 'r') as write_fasta:
	with open(infile+'.strand_switch.faa', 'w') as switch:
		for definition, sequences in SimpleFastaParser(write_fasta):
			if format == "nucl":
				if str(definition).split(" # ")[0].rsplit("_",1)[0] in correct:
					switch.write('>' + (str(definition).split(" # ")[0]+"_"+str(definition).split(" # ")[3].strip("~")) + '\n' + str(sequences) + '\n')
			if format == "prot":
				if str(str(definition).split(" # ")[0].rsplit("_",1)[0]).replace(" ", "$~&").replace('"','^@%') in correct:
					switch.write('>' + (str(definition).split(" # ")[0]+"_"+str(definition).split(" # ")[3].strip("~")).replace(" ", "$~&").replace('"','^@%') + '\n' + str(sequences) + '\n')

with open(infile+'.strand_switch.faa', 'r') as switch:
	switch_list = []
	complete_list = []
	for name, seq in SimpleFastaParser(switch):
		complete_list.append(name.rsplit("_",2)[0])
		complete_list.append(name.rsplit("_",2)[1])
		complete_list.append(name.rsplit("_",2)[2])
	complete_list.append("placeholder")
	n = 0
	genes = 0
	strands = 0
	strand_database = {}
	while n < len(complete_list)-3:
		if complete_list[n] == complete_list[n+3]:
			genes += 1
			if complete_list[n+2] != complete_list[n+5]:
				strands += 1
			n += 3
		if complete_list[n] != complete_list[n+3]:
			genes += 1
			strand_database.update({complete_list[n]:[genes,strands]})
			genes = 0
			strands = 0
			n += 3

with open(infile+'.strand_switch.faa', 'r') as switch:
	with open(infile+'.first_pass_low.faa', 'w') as low:
		with open(infile+'.first_pass_mid.faa', 'w') as mid:
			with open(infile+'.first_pass_high.faa', 'w') as high:
				low_switch = []
				mid_switch = []
				high_switch = []
				for name, seq in SimpleFastaParser(switch):
					denominator = int(strand_database[name.rsplit("_",2)[0]][0])
					numerator = int(strand_database[name.rsplit("_",2)[0]][1])
					if numerator/denominator < 0.05:
						low.write(">" + str(name.rsplit("_",1)[0]) + "\n" + str(seq) + "\n")
						low_switch.append(name)
					elif numerator/denominator < 0.35:
						mid.write(">" + str(name.rsplit("_",1)[0]) + "\n" + str(seq) + "\n")
						mid_switch.append(name)
					elif numerator/denominator >= 0.35:
						high.write(">" + str(name) + "\n" + str(seq) + "\n")
						high_switch.append(name)

if len(high_switch) > 0:
	with open(infile+'.first_pass_high.faa', 'r') as high:
		subprocess.call(['hmmsearch', '--tblout', infile+'.vpfam.hmmtbl', '--noali', '-T', '50', '--cpu', cpu, '-o', infile+'_temp.txt', vpfam_hmm, infile+'.first_pass_high.faa'])
		with open(infile+'.vpfam.hmmtbl', 'r') as vpfam_infile:
			with open(infile+'.vpfam.hmmtbl.temp.txt', 'w') as vpfam_outfile:
				vpfam_outfile.write("protein" + "\t" + "id" + "\t" + "evalue" + "\t" + "score" + '\n')
				for line in vpfam_infile:
					if not line.startswith('#'):
						line = line.split(' ')
						parse = list(filter(None, line))
						vpfam_outfile.write(str(parse[0]) + '\t' + str(parse[3]) + '\t' + str(parse[4]) + '\t' + str(parse[5]) + '\n')
		subprocess.call(['rm', infile+'_temp.txt'])

		with open(infile+'.vpfam.hmmtbl.temp.txt', 'r') as vpfam_sort:
			with open(infile+'.vpfam.hmmtbl.parse.txt', 'w') as vpfam_table:
				table = pd.read_csv(vpfam_sort, sep="\t")
				sort = table.sort_values(by='evalue', ascending=True)
				drop = sort.drop_duplicates(subset='protein', keep='first')
				write = drop.to_csv(vpfam_table, index=False, sep="\t")
		vpfam = infile+".vpfam.hmmtbl.parse.txt"
		subprocess.call(['rm', infile+'.vpfam.hmmtbl.temp.txt'])

##########################  Removing Bacteria  #################################
	with open(vpfam, 'r') as vpfam:
		n = 4
		vpfam = vpfam.read().replace('\n','\t').split('\t')
		if vpfam[-1] == '':
			del vpfam[-1]
		vpfam_list = []
		while n < len(vpfam):
			vpfam_list.append(str(vpfam[n].rsplit("_",2)[0]))
			n += 4
		vpfam_counts = list(set(vpfam_list))
		vpfam_virus = []
		for genome in vpfam_counts:
			num_genes = vpfam_list.count(genome)
			if num_genes > 0:
				vpfam_virus.append(genome)

########################### Generate New Database  #############################
	with open(infile+'.first_pass_high.faa', 'r') as first_pass:
		with open(infile+'.appended_high.faa', 'w') as appended:
			counter = 0
			for name, seq in SimpleFastaParser(first_pass):
				if str(name.rsplit("_",2)[0]) in vpfam_virus:
					appended.write(">" + str(name.rsplit("_",1)[0]) + "\n" + str(seq) + "\n")
					counter += 1

	if len(mid_switch) > 0 and counter > 0:
		cat_file = str('cat ' + infile+'.appended_high.faa' + ' ' + infile+'.first_pass_mid.faa' + ' > ' + infile+'.first_pass.faa')
		subprocess.run(cat_file, shell=True)
		no_contigs = False
	elif counter > 0:
		subprocess.call(['mv', infile+'.appended_high.faa', infile+'.first_pass.faa'])
		no_contigs = False
	else:
		no_contigs = True
#######
elif len(mid_switch) > 0:
	subprocess.call(['mv', infile+'.first_pass_mid.faa', infile+'.first_pass.faa'])
	no_contigs = False
else:
	no_contigs = True

if no_contigs == False:
############################### Run/Parse plasmid hmmsearch  #######################
	subprocess.call(['hmmsearch', '--tblout', infile+'.plasmid.hmmtbl', '--noali', '-T', '50', '--cpu', cpu, '-o', infile+'_temp.txt', plasmid_hmm, infile+'.first_pass.faa'])
	with open(infile+'.plasmid.hmmtbl', 'r') as plasmid_infile:
		with open(infile+'.plasmid.hmmtbl.temp.txt', 'w') as plasmid_outfile:
			plasmid_outfile.write("protein" + "\t" + "id" + "\t" + "evalue" + "\t" + "score" + '\n')
			for line in plasmid_infile:
				if not line.startswith('#'):
					line = line.split(' ')
					parse = list(filter(None, line))
					plasmid_outfile.write(str(parse[0]) + '\t' + str(parse[3]) + '\t' + str(parse[4]) + '\t' + str(parse[5]) + '\n')
	subprocess.call(['rm', infile+'_temp.txt'])

	with open(infile+'.plasmid.hmmtbl.temp.txt', 'r') as plasmid_sort:
		with open(infile+'.plasmid.hmmtbl.parse.txt', 'w') as plasmid_table:
			table = pd.read_csv(plasmid_sort, sep="\t")
			sort = table.sort_values(by='evalue', ascending=True)
			drop = sort.drop_duplicates(subset='protein', keep='first')
			write = drop.to_csv(plasmid_table, index=False, sep="\t")
	plasmid = infile+".plasmid.hmmtbl.parse.txt"
	subprocess.call(['rm', infile+'.plasmid.hmmtbl.temp.txt'])

	##########################  Removing Bacteria  #################################
	with open(plasmid, 'r') as plasmid:
		n = 4
		plasmid = plasmid.read().replace('\n','\t').split('\t')
		if plasmid[-1] == '':
			del plasmid[-1]
		plasmid_list = []
		while n < len(plasmid):
			plasmid_list.append(str(plasmid[n].rsplit("_",2)[0]))
			n += 4
		plasmid_counts = list(set(plasmid_list))
		bacteria = []
		for genome in plasmid_counts:
			num_genes = plasmid_list.count(genome)
			if num_genes >= 3:
				bacteria.append(genome)

	########################### Generate New Database  #############################
	with open(infile+'.first_pass.faa', 'r') as first_pass:
		with open(infile+'.first_pass_mid-high.faa', 'w') as appended:
			counter = 0
			for name, seq in SimpleFastaParser(first_pass):
				if str(name.rsplit("_",1)[0]) not in bacteria:
					appended.write(">" + str(name) + "\n" + str(seq) + "\n")
					counter += 1

	if counter > 0 and len(low_switch) > 0:
		cat_file = str('cat ' + infile+'.first_pass_low.faa' + ' ' + infile+'.first_pass_mid-high.faa' + ' > ' + infile+'.appended.faa')
		subprocess.run(cat_file, shell=True)
	elif counter == 0 and len(low_switch) > 0:
		cat_file = str('mv ' + infile+'.first_pass_low.faa' + ' ' + infile+'.appended.faa')
		subprocess.run(cat_file, shell=True)
	elif counter > 0:
		cat_file = str('mv ' + infile+'.first_pass_mid-high.faa' + ' ' + infile+'.appended.faa')
		subprocess.run(cat_file, shell=True)
	else:
		exit()

else:
	if len(low_switch) == 0:
		exit()
	else:
		subprocess.call(['mv', infile+'.first_pass_low.faa', infile+'.appended.faa'])

rm_first_pass = 'rm '+infile+'.first_pass*.faa'
subprocess.run(rm_first_pass, shell=True)
subprocess.call(['rm', infile+'.strand_switch.faa'])

########################### Generate Master File  ##############################
with open(infile+'.appended.faa', 'r') as write_fasta:
	with open(infile+'.master.txt', 'w') as master:
		master.write("protein" + "\t" + "genome" + "\n")
		counter = 0
		for name, seq in SimpleFastaParser(write_fasta):
			master.write(name + '\t' + str(name.rsplit("_",1)[0]) + '\n')
			counter += 1
if counter == 0:
	exit()

###########################  Run/Parse KEGG hmmsearch ##########################
subprocess.call(['hmmsearch', '--tblout', infile+'.KEGG.hmmtbl', '--noali', '-T', '40', '--cpu', cpu, '-o', infile+'_temp.txt', kegg_hmm, infile+'.appended.faa'])
with open(infile+'.KEGG.hmmtbl', 'r') as kegg_infile:
	with open(infile+'.KEGG.hmmtbl.temp.txt', 'w') as kegg_outfile:
		kegg_outfile.write("protein" + "\t" + "id" + "\t" + "evalue" + "\t" + "score" + '\n')
		for line in kegg_infile:
			if not line.startswith('#'):
				line = line.split(' ')
				parse = list(filter(None, line))
				kegg_outfile.write(str(parse[0]) + '\t' + str(parse[2]) + '\t' + str(parse[4]) + '\t' + str(parse[5]) + '\n')
subprocess.call(['rm', infile+'_temp.txt'])

with open(infile+'.KEGG.hmmtbl.temp.txt', 'r') as kegg_sort:
	with open(infile+'.KEGG.hmmtbl.parse.txt', 'w') as kegg_table:
		table = pd.read_csv(kegg_sort, sep="\t")
		sort = table.sort_values(by='evalue', ascending=True)
		drop = sort.drop_duplicates(subset='protein', keep='first')
		write = drop.to_csv(kegg_table, index=False, sep="\t")
kegg = infile+".KEGG.hmmtbl.parse.txt"
subprocess.call(['rm', infile+'.KEGG.hmmtbl.temp.txt'])

############################### Prophage Search  ###############################

with open(infile+'.master.txt', 'r') as master:
	master = master.read().replace('\n','\t').split('\t')
	with open(str(path)+'temp_kegg_annotations.' + str(base) + '.txt', 'w') as output:
		output.write('protein' + '\t' + 'genome' + '\t' + 'KO' + '\t' + 'KO evalue' + '\t' + 'KO score' + '\t' + 'KO category'  + '\n')

##########################  Creat Dictionaries  ################################
#######  KEGG
		with open(str(kegg), 'r') as kegg:
			n = 4
			kegg = kegg.read().replace('\n','\t').split('\t')
			if kegg[-1] == '':
				del kegg[-1]
			kegg_dict = {}

			while n < len(kegg):
				protein = kegg[n]
				id = kegg[n+1]
				evalue = kegg[n+2]
				score = kegg[n+3]
				combo = str(id + '\t' + evalue + '\t' + score)
				kegg_dict.update({protein:combo})
				n += 4

####### Categories
		with open(categories, 'r', encoding='utf-8-sig') as cat:
			n = 2
			cat = cat.read().replace('\n','\t').split('\t')
			cat_dict = {}

			while n < len(cat):
				id = cat[n]
				value = str((float(cat[n+1])/100))
				cat_dict.update({id:value})
				n += 2

		i = 2
		while i < len(master)-1:
			output.write(master[i]+'\t'+master[i+1])

		##### KEGG
			if master[i] in kegg_dict:
				result = kegg_dict[master[i]]
				result = result.split('\t')
				if result[0] in cat_dict:
					category = cat_dict[result[0]]
				if result[0] not in cat_dict:
					category = "0"
				output.write('\t' + result[0] + '\t' + result[1] + '\t' + result[2] + '\t' + category + '\n')

			if master[i] not in kegg_dict:
				output.write('\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n')
		##### next
			i += 2

with open(str(path)+'temp_kegg_annotations.' + str(base) + '.txt', 'r') as annotations:
####### Create unique genomes List
	annotations = annotations.read().replace('\n','\t').split('\t')
	annotations.append('')
	n = 7
	i = 0
	full_genomes = list()
	genomes = list()
	while n < len(annotations):
		full_genomes.append(annotations[n])
		n += 6
	for item in full_genomes:
		if item not in genomes:
			genomes.append(item)
####### Count annotations per genome
	n = 7
	check1 = 0
	prophage_list = []
	prophage_dict = {}

	for genome in genomes[:-1]:
		check2 = 0
		check1 += 1
		if check1 > 1:
			n -= 6
		total_genes = full_genomes.count(genome)
		gene_number = 1
		genes = []
		frag = 1
		while gene_number <= total_genes:
			genes.append(str(genome)+"_"+str(gene_number))
			gene_number += 1
####### Number of annotations
		counter = 0
		while counter == 0:
			cut = False
			bact_kegg = 0
			blank_kegg = 0
			total_kegg = 0
			zero_kegg = 0

			if n+30 < len(annotations):
				if genome == annotations[n] and genome == annotations[n+18]:
					if annotations[n+4] == '':
						bact_kegg += 1
						blank_kegg += 1
					elif float(annotations[n+4]) <= 0.02:
						bact_kegg += 1
						total_kegg += float(annotations[n+4])
						if float(annotations[n+4]) == 0.0:
							zero_kegg += 1
					elif float(annotations[n+4]) > 0.02:
						bact_kegg -= 1
					if annotations[n+10] == '':
						bact_kegg += 1
						blank_kegg += 1
					elif float(annotations[n+10]) <= 0.02:
						bact_kegg += 1
						total_kegg += float(annotations[n+10])
						if float(annotations[n+10]) == 0.0:
							zero_kegg += 1
					elif float(annotations[n+10]) > 0.02:
						bact_kegg -= 1
					if annotations[n+16] == '':
						bact_kegg += 1
						blank_kegg += 1
					elif float(annotations[n+16]) <= 0.02:
						bact_kegg += 1
						total_kegg += float(annotations[n+16])
						if float(annotations[n+16]) == 0.0:
							zero_kegg += 1
					elif float(annotations[n+16]) > 0.02:
						bact_kegg -= 1
					if annotations[n+22] == '':
						bact_kegg += 1
						blank_kegg += 1
					elif float(annotations[n+22]) <= 0.02:
						bact_kegg += 1
						total_kegg += float(annotations[n+22])
						if float(annotations[n+22]) == 0.0:
							zero_kegg += 1
					elif float(annotations[n+22]) > 0.02:
						bact_kegg -= 1

					if bact_kegg == 4 and blank_kegg <= 1 and total_kegg < 0.06:
						if annotations[n-1] in genes:
							if check2 == 0:
								split2 = int(annotations[n-1].rsplit("_",1)[1])
								fragment = genes[:split2]
							elif check2 > 0:
								split2 = int(annotations[n-7].rsplit("_",1)[1])
								fragment = genes[split1:split2]
							if len(fragment) >= 8:
								prophage_dict.update({str(genome)+"_fragment_"+str(frag):fragment})
								frag += 1
							prophage_list.append(genome)
							prophage_list.append(frag)
							check2 += 1
							split1 = int(annotations[n+17].rsplit("_",1)[1])
							n += 6
							cut = True
					elif zero_kegg == 3:
						if annotations[n-1] in genes:
							if check2 == 0:
								split2 = int(annotations[n-1].rsplit("_",1)[1])
								fragment = genes[:split2]
							elif check2 > 0:
								split2 = int(annotations[n-7].rsplit("_",1)[1])
								fragment = genes[split1:split2]
							if len(fragment) >= 8:
								prophage_dict.update({str(genome)+"_fragment_"+str(frag):fragment})
								frag += 1
							prophage_list.append(genome)
							prophage_list.append(frag)
							check2 += 1
							split1 = int(annotations[n+11].rsplit("_",1)[1])
							n += 6
							cut = True

					elif n+66 < len(annotations):
						if genome == annotations[n] and genome == annotations[n+48]:
							if annotations[n+4] != '' and annotations[n+10] != '' and annotations[n+16] != '' and annotations[n+22] != '' and annotations[n+28] != '' and annotations[n+34] != '' and annotations[n+40] != '' and annotations[n+46] != '' and annotations[n+52] != '':
								if annotations[n-1] in genes:
									if check2 == 0:
										split2 = int(annotations[n-1].rsplit("_",1)[1])
										fragment = genes[:split2]
									elif check2 > 0:
										split2 = int(annotations[n-7].rsplit("_",1)[1])
										fragment = genes[split1:split2]
									if len(fragment) >= 8:
										prophage_dict.update({str(genome)+"_fragment_"+str(frag):fragment})
										frag += 1
									prophage_list.append(genome)
									prophage_list.append(frag)
									check2 += 1
									split1 = int(annotations[n+47].rsplit("_",1)[1])
									n += 6
									cut = True

			elif n+8 < len(annotations):
				if genome == annotations[n] and genome == annotations[n+6]:
					if annotations[n+1] != '' and annotations[n+7] != '':
						if annotations[n+4] == '0' and annotations[n+10] == '0':
							if annotations[n-1] in genes:
								split2 = int(annotations[n-1].rsplit("_",1)[1])
								if check2 == 0:
									fragment = genes[:split2]
								if check2 > 0:
									fragment = genes[split1:split2]
								if len(fragment) >= 8:
									prophage_dict.update({str(genome)+"_fragment_"+str(frag):fragment})
									frag += 1
								prophage_list.append(genome)
								prophage_list.append(frag)
								check2 += 1
								split1 = int(annotations[n-1].rsplit("_",1)[1])
								n += 6
								cut = True
####### next
			if genome != annotations[n]:
				if check2 > 1:
					fragment = genes[split1:]
					if len(fragment) >= 8:
						prophage_dict.update({str(genome)+"_fragment_"+str(frag):fragment})
						frag += 1
					prophage_list.append(genome)
					prophage_list.append(frag)
				counter += 1

			if cut == False:
				n += 6

####### write out new protein file
removed = []
if len(prophage_dict) > 0:
	with open(infile+'.appended.faa', 'r') as read_fasta:
		with open(infile+'.prophage_analysis.faa', 'w') as write_fasta:
			db_dict = {}
			for name, seq in SimpleFastaParser(read_fasta):
				db_dict.update({name:seq})
			for key in prophage_dict:
				n = 0
				while n < len(prophage_dict[key]):
					if prophage_dict[key][n] in db_dict:
						write_fasta.write(">" + str(key) + "_" + str(prophage_dict[key][n].rsplit("_",1)[1]) + "\n" + str(db_dict[prophage_dict[key][n]]) + "\n")
						removed.append(prophage_dict[key][n].rsplit("_",1)[0])
						n += 1
					else:
						n += 1
	removed = list(set(removed))
	with open(infile+'.appended.faa', 'r') as read_fasta:
		with open(infile+'.intact_analysis.faa', 'w') as write_fasta:
			for name, seq in SimpleFastaParser(read_fasta):
				if name.rsplit("_",1)[0] not in removed:
					write_fasta.write(">" + str(name) + "\n" + str(seq) + "\n")

	cat_file = str("cat " + str(infile)+'.intact_analysis.faa' + " " + str(infile)+'.prophage_analysis.faa' + " > " + str(infile)+'.analysis.faa')
	subprocess.run(cat_file, shell=True)
	subprocess.call(['rm', infile+'.intact_analysis.faa', infile+'.prophage_analysis.faa'])

else:
	subprocess.call(['mv', infile+'.appended.faa', infile+'.analysis.faa'])

########################### New Master File ####################################
with open(infile+'.analysis.faa', 'r') as write_fasta:
	with open(infile+'.master.txt', 'w') as master:
		master.write("protein" + "\t" + "genome" + "\n")
		for name, seq in SimpleFastaParser(write_fasta):
			master.write(name + '\t' + str(name.rsplit("_",1)[0]) + '\n')

##########################  Creat Dictionaries  ################################
with open(infile+'.master.txt', 'r') as master:
	master = master.read().replace('\n','\t').split('\t')
	with open(str(path)+'temp_kegg_annotations.' + str(base) + '.txt', 'w') as output:
		output.write('protein' + '\t' + 'genome' + '\t' + 'KO' + '\t' + 'KO evalue' + '\t' + 'KO score' + '\t' + 'KO category'  + '\n')

#######  KEGG
		kegg = infile+".KEGG.hmmtbl.parse.txt"
		with open(kegg, 'r') as kegg:
			n = 4
			kegg = kegg.read().replace('\n','\t').split('\t')
			if kegg[-1] == '':
				del kegg[-1]
			kegg_dict = {}

			while n < len(kegg):
				protein = kegg[n]
				id = kegg[n+1]
				evalue = kegg[n+2]
				score = kegg[n+3]
				combo = str(id + '\t' + evalue + '\t' + score)
				kegg_dict.update({protein:combo})
				n += 4
		i = 2
		while i < len(master)-1:
			fragment_master = str(master[i].rsplit("_",3)[0] + "_" + master[i].rsplit("_",1)[1])
			output.write(master[i]+'\t'+master[i+1])

		##### KEGG
			if "fragment" in str(master[i]):
				if str(fragment_master) in kegg_dict:
					result = kegg_dict[fragment_master]
					result = result.split('\t')
					if result[0] in cat_dict:
						category = cat_dict[result[0]]
					if result[0] not in cat_dict:
						category = "0"
					output.write('\t' + result[0] + '\t' + result[1] + '\t' + result[2] + '\t' + category + '\n')
				if str(fragment_master) not in kegg_dict:
					output.write('\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n')
			##### next
				i += 2

			else:
				if master[i] in kegg_dict:
					result = kegg_dict[master[i]]
					result = result.split('\t')
					if result[0] in cat_dict:
						category = cat_dict[result[0]]
					if result[0] not in cat_dict:
						category = "0"
					output.write('\t' + result[0] + '\t' + result[1] + '\t' + result[2] + '\t' + category + '\n')
				if master[i] not in kegg_dict:
					output.write('\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n')
			##### next
				i += 2

###########################  Calculate Counts  #################################
with open(str(path)+'temp_kegg_annotations.' + str(base) + '.txt', 'r') as annotations:

####### Create unique genomes List
	annotations = annotations.read().replace('\n','\t').split('\t')
	annotations.append('')
	kegg_int_rep_list = ['K07729', 'K03733', 'K04763', 'K14059', 'K21039', 'K01356', 'K18918', 'K07741', 'K21528', 'K06400', 'K01152', 'K07481', 'K07482', 'K07483', 'K07484', 'K07485', 'K07486', 'K07487', 'K07488', 'K07489', 'K07491', 'K07492', 'K07493', 'K07494', 'K07495', 'K07496', 'K07497', 'K07498', 'K07499', 'K18320', 'K23209']
	n = 7
	i = 0
	full_genomes = list()
	genomes = list()
	keep = []
	while n < len(annotations):
		full_genomes.append(annotations[n])
		n += 6
	for item in full_genomes:
		if item not in genomes:
			genomes.append(item)

	####### Count annotations per genome
	n = 7
	check1 = 0
	for genome in genomes[:-1]:
		kegg_int_rep_count = 0
		zero = 0
		check1 += 1
		if check1 > 1:
			n -= 6
		all_kegg = 0
		cat_kegg = 0
		kegg_temp = 0
		kegg_in_row = 0
		total_genes = full_genomes.count(genome)

	####### Number of annotations
		counter = 0
		while counter == 0:
			if n == len(annotations)-1:
				break
			if genome == annotations[n]:
				if annotations[n+1] != '':
					kegg_temp += 1
					all_kegg += 1
					cat_kegg += float(annotations[n+4])
					if annotations[n+4] == "0":
						zero += 1
					if annotations[n+1] in kegg_int_rep_list:
						kegg_int_rep_count += 1
				else:
					kegg_temp = 0
			if kegg_temp >= 9 and zero >= 1:
				kegg_in_row += 1
####### next
			if genome != annotations[n]:
				counter += 1
			n += 6

####### write out protein for genomes that don't have too many KEGGs
		if total_genes >= int(orf_low):
			remove_check = 0
			if all_kegg/total_genes < 0.7 and zero < 15 and total_genes < 15:
				keep.append(genome)
				remove_check += 1
			elif all_kegg/total_genes < 0.5 and zero < 15 and total_genes >= 15:
				keep.append(genome)
				remove_check += 1
			elif all_kegg/total_genes < 0.8 and cat_kegg/all_kegg >= 0.2:
				keep.append(genome)
				remove_check += 1
			if "_fragment_" in str(genome) and total_genes < 8 and remove_check > 0:
				keep.remove(genome)
				remove_check = 0
			if all_kegg >= 8 and remove_check > 0:
				if zero / all_kegg > 0.6 and all_kegg/total_genes > 0.2:
					keep.remove(genome)
					remove_check = 0
			if all_kegg < 8 and all_kegg/total_genes > 0.2 and remove_check > 0:
				if zero / all_kegg > 0.6:
					keep.remove(genome)
					remove_check = 0
			if all_kegg < 8 and all_kegg > 2 and "_fragment_" in str(genome) and remove_check > 0:
				if zero / all_kegg >= 0.6:
					keep.remove(genome)
					remove_check = 0
			if all_kegg > 0 and remove_check > 0:
				if all_kegg/total_genes >= 0.5 and cat_kegg/all_kegg <= 0.2:
					keep.remove(genome)
					remove_check = 0
			if kegg_in_row > 0 and remove_check > 0:
				keep.remove(genome)
				remove_check = 0
			if kegg_int_rep_count > 0 and remove_check == 0:
				if zero/all_kegg < 0.35:
					keep.append(genome)
					remove_check += 1

if len(keep) == 0:
	subprocess.call(['rm', str(path)+'temp_kegg_annotations.' + str(base) + '.txt'])
	exit()

with open(infile+'.analysis.faa', 'r') as read_fasta:
	with open(infile+'.appended.faa', 'w') as write_fasta:
		for name, seq in SimpleFastaParser(read_fasta):
			if str(name).rsplit("_",1)[0] in keep:
				write_fasta.write(">" + str(name) + "\n" + str(seq) + "\n")

subprocess.call(['rm', infile+'.analysis.faa'])
subprocess.call(['rm', str(path)+'temp_kegg_annotations.' + str(base) + '.txt'])

########################### Generate Master File  ##############################
with open(infile+'.appended.faa', 'r') as write_fasta:
	with open(infile+'.master.txt', 'w') as master:
		master.write("protein" + "\t" + "genome" + "\n")
		for name, seq in SimpleFastaParser(write_fasta):
			master.write(name + '\t' + str(name.rsplit("_",1)[0]) + '\n')

###########################  Run/Parse Pfam hmmsearch ##########################
subprocess.call(['hmmsearch', '--tblout', infile+'.Pfam.hmmtbl', '--noali', '-T', '40', '--cpu', cpu, '-o', infile+'_temp.txt', pfam_hmm, infile+'.appended.faa'])
with open(infile+'.Pfam.hmmtbl', 'r') as pfam_infile:
	with open(infile+'.Pfam.hmmtbl.temp.txt', 'w') as pfam_outfile:
		pfam_outfile.write("protein" + "\t" + "id" + "\t" + "evalue" + "\t" + "score" + '\n')
		for line in pfam_infile:
			if not line.startswith('#'):
				line = line.split(' ')
				parse = list(filter(None, line))
				pfam_outfile.write(str(parse[0]) + '\t' + str(parse[3]) + '\t' + str(parse[4]) + '\t' + str(parse[5]) + '\n')
subprocess.call(['rm', infile+'_temp.txt'])

with open(infile+'.Pfam.hmmtbl.temp.txt', 'r') as pfam_sort:
	with open(infile+'.Pfam.hmmtbl.parse.txt', 'w') as pfam_table:
		table = pd.read_csv(pfam_sort, sep="\t")
		sort = table.sort_values(by='evalue', ascending=True)
		drop = sort.drop_duplicates(subset='protein', keep='first')
		write = drop.to_csv(pfam_table, index=False, sep="\t")
pfam = infile+".Pfam.hmmtbl.parse.txt"
subprocess.call(['rm', infile+'.Pfam.hmmtbl.temp.txt'])


with open(infile+'.master.txt', 'r') as master:
	master = master.read().replace('\n','\t').split('\t')
	with open(str(path)+'temp_pfam_annotations.' + str(base) + '.txt', 'w') as output:
		output.write('protein' + '\t' + 'genome' + '\t' + 'Pfam' + '\t' + 'Pfam evalue' + '\t' + 'Pfam score' + '\t' + 'Pfam category'  + '\n')

#######  pfam
		pfam = infile+".Pfam.hmmtbl.parse.txt"
		with open(pfam, 'r') as pfam:
			n = 4
			pfam = pfam.read().replace('\n','\t').split('\t')
			if pfam[-1] == '':
				del pfam[-1]
			pfam_dict = {}
			while n < len(pfam):
				protein = pfam[n]
				id = pfam[n+1]
				evalue = pfam[n+2]
				score = pfam[n+3]
				combo = str(id + '\t' + evalue + '\t' + score)
				pfam_dict.update({protein:combo})
				n += 4
		i = 2
		while i < len(master)-1:
			output.write(master[i]+'\t'+master[i+1])

			if master[i] in pfam_dict:
				result = pfam_dict[master[i]]
				result = result.split('\t')
				if result[0] in cat_dict:
					category = cat_dict[result[0]]
				if result[0] not in cat_dict:
					category = "0"
				output.write('\t' + result[0] + '\t' + result[1] + '\t' + result[2] + '\t' + category + '\n')
			if master[i] not in pfam_dict:
				output.write('\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n')
		##### next
			i += 2

###########################  Calculate Counts  #################################
with open(str(path)+'temp_pfam_annotations.' + str(base) + '.txt', 'r') as annotations:
####### Create unique genomes List
	annotations = annotations.read().replace('\n','\t').split('\t')
	annotations.append('')
	pfam_int_rep_list = ['PF00239.21', 'PF00589.22', 'PF00665.26', 'PF00872.18', 'PF01385.19', 'PF01526.17', 'PF01527.20', 'PF01548.17', 'PF01609.21', 'PF01610.17', 'PF01726.16', 'PF01797.16', 'PF01818.17', 'PF02061.16', 'PF02316.16', 'PF02899.17', 'PF02914.15', 'PF03050.14', 'PF03374.14', 'PF03400.13', 'PF04754.12', 'PF04761.12', 'PF05269.11', 'PF05598.11', 'PF05973.14', 'PF06543.12', 'PF07022.13', 'PF07282.11', 'PF07508.13', 'PF09003.10', 'PF09299.11', 'PF09588.10', 'PF09669.10', 'PF10551.9', 'PF12167.8', 'PF12472.8', 'PF12760.7', 'PF13022.6', 'PF13333.6', 'PF13340.6', 'PF13356.6', 'PF13408.6', 'PF13586.6', 'PF13701.6', 'PF13737.6', 'PF14657.6', 'PF14659.6', 'PF16452.5', 'PF16795.5']
	n = 7
	i = 0
	full_genomes = list()
	genomes = list()
	keep = []
	while n < len(annotations):
		full_genomes.append(annotations[n])
		n += 6
	for item in full_genomes:
		if item not in genomes:
			genomes.append(item)

	####### Count annotations per genome
	n = 7
	check1 = 0
	for genome in genomes[:-1]:
		pfam_int_rep_count = 0
		zero = 0
		check1 += 1
		if check1 > 1:
			n -= 6
		all_pfam = 0
		cat_pfam = 0
		total_genes = full_genomes.count(genome)

	####### Number of annotations
		counter = 0
		while counter == 0:
			if n == len(annotations)-1:
				break
			if genome == annotations[n]:
				if annotations[n+1] != '':
					all_pfam += 1
					cat_pfam += float(annotations[n+4])
					if annotations[n+4] == "0":
						zero += 1
					if annotations[n+1] in pfam_int_rep_list:
						pfam_int_rep_count += 1
####### next
			if genome != annotations[n]:
				counter += 1
			n += 6
####### write out protein for genomes that don't have too many pfams
		if total_genes >= int(orf_low):
			remove_check = 0
			if all_pfam <= 15:
				keep.append(genome)
				remove_check += 1
			elif all_pfam/total_genes <= 0.6:
				keep.append(genome)
				remove_check += 1
			elif all_pfam/total_genes > 0.6 and cat_pfam/all_pfam >= 0.15:
				keep.append(genome)
				remove_check += 1
			if zero >= 15 and remove_check > 0:
				keep.remove(genome)
				remove_check = 0
			if all_pfam > 8 and remove_check > 0:
				if zero / all_pfam > 0.5 and all_pfam/total_genes > 0.2:
					keep.remove(genome)
					remove_check = 0
			if all_pfam > 0 and remove_check > 0:
				if zero / all_pfam > 0.6 and all_pfam/total_genes > 0.2:
					keep.remove(genome)
					remove_check = 0
			if pfam_int_rep_count > 0 and "_fragment_" in str(genome) and remove_check == 0:
				if zero/all_pfam < 0.35 and all_pfam > 4:
					keep.append(genome)
					remove_check += 1
				elif zero/all_pfam <= 0.5 and all_pfam <= 5 and total_genes > 10:
					keep.append(genome)
					remove_check += 1

if len(keep) == 0:
	subprocess.call(['rm', str(path)+'temp_pfam_annotations.' + str(base) + '.txt'])
	exit()

########################### Generate New Database  #############################
with open(infile+'.appended.faa', 'r') as read_fasta:
	with open(infile+'.pass.faa', 'w') as write_fasta:
		for name, seq in SimpleFastaParser(read_fasta):
			if str(name).rsplit("_",1)[0] in keep:
				write_fasta.write(">" + str(name) + "\n" + str(seq) + "\n")

subprocess.call(['rm', str(path)+'temp_pfam_annotations.' + str(base) + '.txt'])

########################### Generate Master File  ##############################
with open(infile+'.pass.faa', 'r') as write_fasta:
	with open(infile+'.master.txt', 'w') as master:
		master.write("protein" + "\t" + "genome" + "\n")
		for name, seq in SimpleFastaParser(write_fasta):
			master.write(name + '\t' + str(name.rsplit("_",1)[0]) + '\n')

###########################  Run/Parse VOG hmmsearch ##########################

subprocess.call(['hmmsearch', '--tblout', infile+'.VOG.hmmtbl', '--noali', '-T', '40', '--cpu', cpu, '-o', infile+'_temp.txt', vog_hmm, infile+'.pass.faa'])
with open(infile+'.VOG.hmmtbl', 'r') as kegg_infile:
	with open(infile+'.VOG.hmmtbl.temp.txt', 'w') as kegg_outfile:
		kegg_outfile.write("protein" + "\t" + "id" + "\t" + "evalue" + "\t" + "score" + '\n')
		for line in kegg_infile:
			if not line.startswith('#'):
				line = line.split(' ')
				parse = list(filter(None, line))
				kegg_outfile.write(str(parse[0]) + '\t' + str(parse[2]) + '\t' + str(parse[4]) + '\t' + str(parse[5]) + '\n')
subprocess.call(['rm', infile+'_temp.txt'])

with open(infile+'.VOG.hmmtbl.temp.txt', 'r') as kegg_sort:
	with open(infile+'.VOG.hmmtbl.parse.txt', 'w') as kegg_table:
		table = pd.read_csv(kegg_sort, sep="\t")
		sort = table.sort_values(by='evalue', ascending=True)
		drop = sort.drop_duplicates(subset='protein', keep='first')
		write = drop.to_csv(kegg_table, index=False, sep="\t")
vog = infile+".VOG.hmmtbl.parse.txt"
subprocess.call(['rm', infile+'.VOG.hmmtbl.temp.txt'])

############################### Open Files  ####################################
with open(infile+'.master.txt', 'r') as master:
	master = master.read().replace('\n','\t').split('\t')
	with open(str(path)+'temp1_VIBRANT_annotations.' + str(base) + '.txt', 'w') as output:

##########################  Creat Dictionaries  ################################
#######  VOG
		with open(vog, 'r') as vog:
			n = 4
			vog = vog.read().replace('\n','\t').split('\t')
			if vog[-1] == '':
				del vog[-1]
			vog_dict = {}
			while n < len(vog):
				protein = vog[n]
				id = vog[n+1]
				evalue = vog[n+2]
				score = vog[n+3]
				combo = str(id + '\t' + evalue + '\t' + score)
				vog_dict.update({protein:combo})
				n += 4

###########################  Write Annotations  ################################
		i = 2
		while i < len(master)-1:
			fragment_master = str(master[i].rsplit("_",3)[0] + "_" + master[i].rsplit("_",1)[1])
			output.write(master[i]+'\t'+master[i+1])

##### KEGG
			if "fragment" in str(master[i]):
				if str(fragment_master) in kegg_dict:
					result = kegg_dict[fragment_master]
					result = result.split('\t')
					if result[0] in cat_dict:
						category = cat_dict[result[0]]
					elif result[0] not in cat_dict:
						category = "0"
					output.write('\t' + result[0] + '\t' + result[1] + '\t' + result[2] + '\t' + category)

				elif str(fragment_master) not in kegg_dict:
					output.write('\t' + '' + '\t' + '' + '\t' + '' + '\t' + '')

			elif "fragment" not in str(master[i]):
				if master[i] in kegg_dict:
					result = kegg_dict[master[i]]
					result = result.split('\t')
					if result[0] in cat_dict:
						category = cat_dict[result[0]]
					elif result[0] not in cat_dict:
						category = "0"
					output.write('\t' + result[0] + '\t' + result[1] + '\t' + result[2] + '\t' + category)
				elif master[i] not in kegg_dict:
					output.write('\t' + '' + '\t' + '' + '\t' + '' + '\t' + '')

##### Pfam
			if master[i] in pfam_dict:
				result = pfam_dict[master[i]]
				result = result.split('\t')
				if result[0] in cat_dict:
					category = cat_dict[result[0]]
				elif result[0] not in cat_dict:
					category = "0"
				output.write('\t' + result[0] + '\t' + result[1] + '\t' + result[2] + '\t' + category)
			elif master[i] not in pfam_dict:
				output.write('\t' + '' + '\t' + '' + '\t' + '' + '\t' + '')

##### VOG
			if master[i] in vog_dict:
				result = vog_dict[master[i]]
				result = result.split('\t')
				if result[0] in cat_dict:
					category = cat_dict[result[0]]
				elif result[0] not in cat_dict:
					category = "0"
				output.write('\t' + result[0] + '\t' + result[1] + '\t' + result[2] + '\t' + category + '\n')
			elif master[i] not in vog_dict:
				output.write('\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n')
##### next
			i += 2

###########################  Calculate Counts  #################################
with open(str(path)+'temp1_VIBRANT_annotations.' + str(base) + '.txt', 'r') as annotations:
	with open(str(path)+'temp_VIBRANT_results.' + str(base) + '.txt', 'w') as results:
		with open(str(path)+'unmodified_VIBRANT_results.' + str(base) + '.txt', 'w') as unmod_results:
			results.write('genome' + '\t' + 'total genes'  + '\t' + 'all KEGG'  + '\t' + 'category KEGG'  + '\t' + 'all Pfam'  + '\t' + 'category vPfam'  + '\t' + 'all VOG'  + '\t' + 'category VOG'  + '\t' + 'KEGG int-rep'  + '\t' + 'KEGG zero'  + '\t' + 'Pfam int-rep'  + '\t' + 'Pfam zero'  + '\t' + 'VOG redoxin'  + '\t' + 'VOG rec-tran'  + '\t' + 'VOG int'  + '\t' + 'VOG RnR'  + '\t' + 'VOG DNA'  + '\t' + 'KEGG restriction check' + '\t' + 'KEGG toxin check' + '\t' + 'VOG special' + '\t' + 'annotation check'  + '\t' + 'p_v check'  + '\t' + 'p_k check'  + '\t' + 'k_v check'  + '\t' + 'k check'  + '\t' + 'p check'  + '\t' + 'v check'  + '\t' + 'h check' + '\n')
			counter_check = 0
			temp_dict = {}
			temp_list = []
	####### Create unique genomes List
			annotations = annotations.read().replace('\n','\t').split('\t')
			n = 1
			i = 0
			full_genomes = list()
			genomes = list()
			annotations = annotations[:-1]
			while n < len(annotations):
				full_genomes.append(annotations[n])
				n += 14
			for item in full_genomes:
				if item not in genomes:
					genomes.append(item)
			n = 0
			while n < len(annotations):
				annotations_result = str(annotations[n+1]) + "\t" + str(annotations[n+2]) + "\t" + str(annotations[n+3]) + "\t" + str(annotations[n+4]) + "\t" + str(annotations[n+5]) + "\t" + str(annotations[n+6]) + "\t" + str(annotations[n+7]) + "\t" + str(annotations[n+8]) + "\t" + str(annotations[n+9]) + "\t" + str(annotations[n+10]) + "\t" + str(annotations[n+11]) + "\t" + str(annotations[n+12]) + "\t" + str(annotations[n+13])
				temp_dict.update({str(annotations[n]):annotations_result})
				temp_list.append(str(annotations[n]))
				n += 14
			n = 1
			check = 0
			final_temp = []
			ORFs = []
			virus_pred = []
			final_proteins = []
			machine_check = 0

			redoxin  = ['VOG00036', 'VOG00051', 'VOG02803', 'VOG05234', 'VOG11831', 'VOG11834', 'VOG00102'] #also includes methytransferases
			recombinase_transposase = ['VOG11335', 'VOG04941', 'VOG00520', 'VOG02659', 'VOG09005', 'VOG11308', 'VOG00022', 'VOG00758', 'VOG02220', 'VOG08467', 'VOG11977', 'VOG12822', 'VOG17298', 'VOG17662', 'VOG17783', 'VOG19530', 'VOG20981', 'VOG21677', 'VOG22616', 'VOG00654', 'VOG01045', 'VOG10286', 'VOG22418', 'VOG01110', 'VOG04672', 'VOG11317', 'VOG04769', 'VOG14097', 'VOG22615', 'VOG00008', 'VOG11317', 'VOG00449'] #also includes select transcriptional activators, plasmid and exonuclease proteins
			kegg_int_rep_list = ['K07729', 'K03733', 'K04763', 'K14059', 'K21039', 'K01356', 'K18918', 'K07741', 'K21528', 'K06400', 'K01152', 'K07481', 'K07482', 'K07483', 'K07484', 'K07485', 'K07486', 'K07487', 'K07488', 'K07489', 'K07491', 'K07492', 'K07493', 'K07494', 'K07495', 'K07496', 'K07497', 'K07498', 'K07499', 'K18320', 'K23209']
			pfam_int_rep_list = ['PF00239.21', 'PF00589.22', 'PF00665.26', 'PF00872.18', 'PF01385.19', 'PF01526.17', 'PF01527.20', 'PF01548.17', 'PF01609.21', 'PF01610.17', 'PF01726.16', 'PF01797.16', 'PF01818.17', 'PF02061.16', 'PF02316.16', 'PF02899.17', 'PF02914.15', 'PF03050.14', 'PF03374.14', 'PF03400.13', 'PF04754.12', 'PF04761.12', 'PF05269.11', 'PF05598.11', 'PF05973.14', 'PF06543.12', 'PF07022.13', 'PF07282.11', 'PF07508.13', 'PF09003.10', 'PF09299.11', 'PF09588.10', 'PF09669.10', 'PF10551.9', 'PF12167.8', 'PF12472.8', 'PF12760.7', 'PF13022.6', 'PF13333.6', 'PF13340.6', 'PF13356.6', 'PF13408.6', 'PF13586.6', 'PF13701.6', 'PF13737.6', 'PF14657.6', 'PF14659.6', 'PF16452.5', 'PF16795.5']
			integrase = ['VOG00041', 'VOG15133', 'VOG20969', 'VOG02658', 'VOG04024', 'VOG01778', 'VOG02371']
			RnR = ['VOG00068', 'VOG00120', 'VOG00313', 'VOG00398', 'VOG01019', 'VOG03262', 'VOG09395', 'VOG21278', 'VOG15771', 'VOG00688', 'VOG00715'] # and thymidylate synthase
			DNA = ['VOG02853', 'VOG07917', 'VOG24666', 'VOG03664', 'VOG04691', 'VOG21077', 'VOG00821', 'VOG00163', 'VOG11329', 'VOG11335', 'VOG10269', 'VOG11468', 'VOG00817', 'VOG01029', 'VOG00031', 'VOG01431', 'VOG03850', 'VOG00038', 'VOG00108', 'VOG01695', 'VOG00399', 'VOG00654', 'VOG01045', 'VOG00404', 'VOG00029', 'VOG00056', 'VOG00079', 'VOG00084', 'VOG00109', 'VOG00272', 'VOG00329', 'VOG00356', 'VOG00422', 'VOG00474', 'VOG00476', 'VOG00532', 'VOG00614', 'VOG00674', 'VOG00709', 'VOG00724', 'VOG00728', 'VOG00815', 'VOG00863', 'VOG00895', 'VOG00908', 'VOG00935', 'VOG01042', 'VOG01137', 'VOG01374', 'VOG01384', 'VOG01422', 'VOG01450', 'VOG01519', 'VOG01585', 'VOG01605', 'VOG01637', 'VOG01664', 'VOG01697', 'VOG01763', 'VOG01852', 'VOG01895', 'VOG01901', 'VOG02015', 'VOG02108', 'VOG02209', 'VOG02395', 'VOG02401', 'VOG02569', 'VOG02821', 'VOG02851', 'VOG02874', 'VOG03009', 'VOG03044', 'VOG03307', 'VOG03317', 'VOG03439', 'VOG03479', 'VOG03586', 'VOG03953', 'VOG04002', 'VOG04063', 'VOG04302', 'VOG04652', 'VOG04729', 'VOG04850', 'VOG04962', 'VOG04998', 'VOG05011', 'VOG05034', 'VOG05091', 'VOG05095', 'VOG05100', 'VOG05105', 'VOG05115', 'VOG05119', 'VOG05300', 'VOG05317', 'VOG05323', 'VOG05344', 'VOG05382', 'VOG05406', 'VOG05693', 'VOG05802', 'VOG05966', 'VOG05972', 'VOG06283', 'VOG07014', 'VOG07060', 'VOG07770', 'VOG08033', 'VOG08045', 'VOG08144', 'VOG08273', 'VOG08380', 'VOG08409', 'VOG08440', 'VOG08827', 'VOG09260', 'VOG09702', 'VOG09774', 'VOG10038', 'VOG10089', 'VOG10101', 'VOG10140', 'VOG10207', 'VOG10235', 'VOG10289', 'VOG11125', 'VOG11285', 'VOG11351', 'VOG11729', 'VOG11958', 'VOG12206', 'VOG12269', 'VOG14398', 'VOG14781', 'VOG17046', 'VOG17463', 'VOG17529', 'VOG18148', 'VOG18226', 'VOG18480', 'VOG18481', 'VOG18702', 'VOG18819', 'VOG18988', 'VOG19293', 'VOG19373', 'VOG19571', 'VOG19864', 'VOG20018', 'VOG20102', 'VOG20254', 'VOG20418', 'VOG20985', 'VOG21265', 'VOG21353', 'VOG21722', 'VOG21964', 'VOG22116', 'VOG22861', 'VOG22863', 'VOG23467', 'VOG23468', 'VOG23717', 'VOG24749', 'VOG24750', 'VOG24981', 'VOG00393', 'VOG24404', 'VOG10258', 'VOG00557', 'VOG02455', 'VOG03269', 'VOG06969', 'VOG03619', 'VOG12136', 'VOG19468', 'VOG04560', 'VOG15690', 'VOG00222', 'VOG00157', 'VOG01469', 'VOG02646', 'VOG12013', 'VOG08321', 'VOG12635', 'VOG03573', 'VOG00140', 'VOG16387', 'VOG11378', 'VOG24766', 'VOG10041', 'VOG09778', 'VOG03513', 'VOG09828', 'VOG01256', 'VOG01464', 'VOG02533', 'VOG05391', 'VOG07031', 'VOG00098']
			special = ['VOG17265', 'VOG23378', 'VOG03971', 'VOG19974', 'VOG21494', 'VOG19846', 'VOG21010', 'VOG23255', 'VOG23256', 'VOG08837', 'VOG08834', 'VOG21336', 'VOG19582', 'VOG16059', 'VOG16060', 'VOG03412', 'VOG02335', 'VOG23390', 'VOG24774', 'VOG24789', 'VOG23939', 'VOG19698', 'VOG20700', 'VOG21766', 'VOG24620', 'VOG21782', 'VOG24791', 'VOG23557', 'VOG11305', 'VOG22164', 'VOG02628', 'VOG15962', 'VOG17611', 'VOG18579', 'VOG22447', 'VOG23134', 'VOG10066', 'VOG04780', 'VOG18550', 'VOG18551', 'VOG23079', 'VOG23391', 'VOG24118', 'VOG21998', 'VOG06315', 'VOG10278', 'VOG16253', 'VOG19887', 'VOG21123', 'VOG02491', 'VOG19122', 'VOG03884', 'VOG15965', 'VOG11324', 'VOG22086', 'VOG04393', 'VOG19193', 'VOG20183', 'VOG04197', 'VOG19981', 'VOG14348', 'VOG22444', 'VOG04150', 'VOG18810', 'VOG11983', 'VOG12138', 'VOG10279', 'VOG19484', 'VOG23315', 'VOG23490', 'VOG00528', 'VOG16358', 'VOG04297', 'VOG04356', 'VOG23566', 'VOG22745', 'VOG03721', 'VOG11306', 'VOG19832', 'VOG17540', 'VOG00955', 'VOG04139', 'VOG06441', 'VOG06443', 'VOG07385', 'VOG07388', 'VOG07473', 'VOG08323', 'VOG08324', 'VOG08325', 'VOG08326', 'VOG08327', 'VOG08328', 'VOG08340', 'VOG08343', 'VOG14347', 'VOG15036', 'VOG15928', 'VOG15935', 'VOG15939', 'VOG15963', 'VOG19975', 'VOG21139', 'VOG22650', 'VOG12902', 'VOG01051', 'VOG19855', 'VOG11420', 'VOG15930', 'VOG08329', 'VOG02389', 'VOG01861', 'VOG01052', 'VOG04138', 'VOG15936', 'VOG08336', 'VOG15934', 'VOG08865', 'VOG21079', 'VOG23556', 'VOG08341', 'VOG21085', 'VOG15964', 'VOG23077', 'VOG15933', 'VOG15938', 'VOG08332', 'VOG08335', 'VOG02193', 'VOG16504', 'VOG21081', 'VOG04415', 'VOG09120', 'VOG21080', 'VOG08333', 'VOG21084', 'VOG21086', 'VOG04907', 'VOG08334', 'VOG04822', 'VOG04228', 'VOG21088', 'VOG03844', 'VOG04764', 'VOG21082', 'VOG04935', 'VOG03204', 'VOG01734', 'VOG00296', 'VOG03147', 'VOG06128', 'VOG06289', 'VOG17834', 'VOG23226', 'VOG13979', 'VOG18491', 'VOG24613', 'VOG06316', 'VOG06375', 'VOG06213', 'VOG07922', 'VOG02468', 'VOG21446', 'VOG12267', 'VOG03194', 'VOG01505', 'VOG06374', 'VOG06140', 'VOG15391', 'VOG03520', 'VOG18860', 'VOG19382', 'VOG11975', 'VOG06137', 'VOG15761', 'VOG02145', 'VOG04079', 'VOG07639', 'VOG06285', 'VOG12899', 'VOG23376', 'VOG24353', 'VOG18667', 'VOG00367', 'VOG00227', 'VOG19851', 'VOG23555', 'VOG12225', 'VOG23374', 'VOG01647', 'VOG22306', 'VOG21440', 'VOG22879', 'VOG20648', 'VOG24612', 'VOG01694', 'VOG11430', 'VOG11302', 'VOG03777', 'VOG13935', 'VOG04522', 'VOG10851', 'VOG00992', 'VOG17613', 'VOG02117', 'VOG03059', 'VOG03404', 'VOG03926', 'VOG19265', 'VOG19587', 'VOG23227', 'VOG24068', 'VOG23257', 'VOG01346', 'VOG03399', 'VOG01179', 'VOG01973', 'VOG19429', 'VOG05705', 'VOG09906', 'VOG11707', 'VOG21325', 'VOG10212', 'VOG06218', 'VOG11404', 'VOG00156', 'VOG00293', 'VOG00541', 'VOG00564', 'VOG01022', 'VOG01143', 'VOG01447', 'VOG01481', 'VOG01933', 'VOG01941', 'VOG02050', 'VOG02069', 'VOG02360', 'VOG02603', 'VOG02784', 'VOG03114', 'VOG03170', 'VOG03254', 'VOG03303', 'VOG03534', 'VOG04291', 'VOG12621', 'VOG14566', 'VOG14730', 'VOG16449', 'VOG17449', 'VOG17452', 'VOG18314', 'VOG19081', 'VOG19271', 'VOG19806', 'VOG20055', 'VOG20145', 'VOG20156', 'VOG21749', 'VOG21991', 'VOG22176', 'VOG22825', 'VOG23189', 'VOG23575', 'VOG23576', 'VOG23927', 'VOG24826', 'VOG02754', 'VOG02419', 'VOG01341', 'VOG03084', 'VOG24311', 'VOG11764', 'VOG18132', 'VOG23129', 'VOG16767', 'VOG01035', 'VOG17546', 'VOG00257', 'VOG19860', 'VOG02944', 'VOG17022', 'VOG18956', 'VOG00101', 'VOG04298', 'VOG00651', 'VOG17451', 'VOG12616', 'VOG04398', 'VOG12622', 'VOG01720', 'VOG03069', 'VOG00929', 'VOG01372', 'VOG11405', 'VOG03383', 'VOG18346', 'VOG03470', 'VOG23090', 'VOG11456', 'VOG11458', 'VOG03487', 'VOG03143', 'VOG11457', 'VOG22335', 'VOG18264', 'VOG04815', 'VOG12907', 'VOG12072', 'VOG01537', 'VOG13937', 'VOG01910', 'VOG01795', 'VOG02823', 'VOG11683', 'VOG21902', 'VOG02067', 'VOG19992', 'VOG01241', 'VOG12227', 'VOG09900', 'VOG12228', 'VOG11088', 'VOG11982', 'VOG04898', 'VOG19579', 'VOG09336', 'VOG01975', 'VOG07839', 'VOG15568', 'VOG13642', 'VOG15567', 'VOG21165', 'VOG02813', 'VOG03014', 'VOG03217', 'VOG10214', 'VOG10242', 'VOG17242', 'VOG23908', 'VOG10146', 'VOG00489', 'VOG01633', 'VOG09941', 'VOG11985', 'VOG17915', 'VOG02782', 'VOG21767', 'VOG17917', 'VOG01513', 'VOG21669', 'VOG10251', 'VOG10252', 'VOG12917', 'VOG10172', 'VOG17900', 'VOG06155', 'VOG01389', 'VOG10243', 'VOG11711', 'VOG17926', 'VOG15597', 'VOG17914', 'VOG12918', 'VOG23907', 'VOG02413', 'VOG01350', 'VOG02080', 'VOG01771', 'VOG10174', 'VOG00501', 'VOG18043', 'VOG23031', 'VOG23235', 'VOG23375', 'VOG13980', 'VOG18538', 'VOG21509', 'VOG17756', 'VOG10185', 'VOG13943', 'VOG15967', 'VOG18134', 'VOG10265', 'VOG14735', 'VOG24473', 'VOG15541', 'VOG03699', 'VOG24632', 'VOG14594', 'VOG13649', 'VOG10325', 'VOG10538', 'VOG12146', 'VOG19196', 'VOG19985', 'VOG10280', 'VOG22425', 'VOG19483', 'VOG18036', 'VOG12921', 'VOG16362', 'VOG10264', 'VOG01434', 'VOG13459', 'VOG22145', 'VOG21505', 'VOG02996', 'VOG16532', 'VOG04492', 'VOG20016', 'VOG20963', 'VOG22173', 'VOG24447', 'VOG00549', 'VOG11341', 'VOG11344', 'VOG15577', 'VOG23471', 'VOG23462', 'VOG12900', 'VOG20627', 'VOG10260', 'VOG10188', 'VOG00765', 'VOG10147', 'VOG12432', 'VOG01376', 'VOG01809', 'VOG00090', 'VOG03224', 'VOG02684', 'VOG01873', 'VOG01543', 'VOG23273', 'VOG12190', 'VOG10104', 'VOG03310', 'VOG03354', 'VOG01249', 'VOG11459', 'VOG01806', 'VOG11223', 'VOG10106', 'VOG11579', 'VOG11581', 'VOG02473', 'VOG01629', 'VOG02333', 'VOG11583', 'VOG11585', 'VOG02914', 'VOG08331', 'VOG03803', 'VOG10263', 'VOG10262', 'VOG01105', 'VOG11582', 'VOG10261', 'VOG15579', 'VOG02539', 'VOG01321', 'VOG19840', 'VOG15362', 'VOG04533', 'VOG02548', 'VOG20766', 'VOG03436', 'VOG01591', 'VOG11765', 'VOG00959', 'VOG11303', 'VOG24545', 'VOG24825', 'VOG23935', 'VOG02701', 'VOG02490', 'VOG01437', 'VOG01551', 'VOG11142', 'VOG18035', 'VOG00319', 'VOG01793', 'VOG11125', 'VOG11089', 'VOG00067', 'VOG10266', 'VOG00637', 'VOG00748', 'VOG00988', 'VOG00639', 'VOG00364', 'VOG01829', 'VOG01954', 'VOG00095', 'VOG01944', 'VOG01489', 'VOG12605', 'VOG00678', 'VOG04497', 'VOG00192', 'VOG13487', 'VOG01567', 'VOG02076', 'VOG01483', 'VOG03356', 'VOG00700', 'VOG00720', 'VOG00922', 'VOG00167', 'VOG14056', 'VOG14059', 'VOG14055', 'VOG14054', 'VOG01405', 'VOG14051', 'VOG03000', 'VOG14061', 'VOG00710', 'VOG00547', 'VOG00359', 'VOG02998', 'VOG01571', 'VOG02790', 'VOG23481', 'VOG01003', 'VOG11749', 'VOG19498', 'VOG02752', 'VOG00871', 'VOG01495', 'VOG23339', 'VOG00046', 'VOG00530', 'VOG00177', 'VOG10015', 'VOG09766', 'VOG11785', 'VOG01021', 'VOG01027', 'VOG10954', 'VOG10953', 'VOG00460', 'VOG04756', 'VOG04499', 'VOG00681', 'VOG01398', 'VOG10478', 'VOG10026', 'VOG01666', 'VOG02486', 'VOG00106', 'VOG02984', 'VOG00923', 'VOG02623', 'VOG00392', 'VOG00534', 'VOG02574', 'VOG09790', 'VOG19369', 'VOG09812', 'VOG10752', 'VOG02457', 'VOG09829', 'VOG10761', 'VOG00132', 'VOG00334', 'VOG00186', 'VOG00360', 'VOG02510', 'VOG01904', 'VOG09805', 'VOG09806', 'VOG01488', 'VOG02498', 'VOG00721', 'VOG03347', 'VOG01104', 'VOG01630', 'VOG03107', 'VOG10125', 'VOG01889', 'VOG10754', 'VOG10758', 'VOG02435', 'VOG00419', 'VOG00608', 'VOG10316', 'VOG03023', 'VOG00644', 'VOG10950', 'VOG01039', 'VOG03428', 'VOG00179', 'VOG00305', 'VOG00743', 'VOG01162', 'VOG02913', 'VOG03429', 'VOG01797', 'VOG01359', 'VOG11082', 'VOG00218', 'VOG00149', 'VOG10315', 'VOG03272', 'VOG00951', 'VOG00699', 'VOG11949', 'VOG00283', 'VOG00698', 'VOG11950', 'VOG00410', 'VOG00793', 'VOG00097', 'VOG02938', 'VOG01397', 'VOG01966', 'VOG03252', 'VOG10951', 'VOG00230', 'VOG01803', 'VOG02869', 'VOG02762', 'VOG00912', 'VOG00129', 'VOG10317', 'VOG03064', 'VOG00125', 'VOG00381', 'VOG00687', 'VOG00751', 'VOG03528', 'VOG02524', 'VOG09814', 'VOG00744', 'VOG01880', 'VOG00537', 'VOG00686', 'VOG00991', 'VOG09809', 'VOG00066', 'VOG01144', 'VOG00819', 'VOG10319', 'VOG00158', 'VOG02516', 'VOG00771', 'VOG00295', 'VOG00981', 'VOG02615', 'VOG00069', 'VOG02877', 'VOG00243', 'VOG03165', 'VOG00324', 'VOG10433', 'VOG00241', 'VOG00308', 'VOG00677', 'VOG00326', 'VOG00727', 'VOG00470', 'VOG01125', 'VOG00236', 'VOG00633', 'VOG01646', 'VOG00213', 'VOG00064', 'VOG09810', 'VOG01297', 'VOG00047', 'VOG00921', 'VOG00321', 'VOG01965', 'VOG01802', 'VOG04946']
			toxin = ['K21498', 'K21497', 'K21496', 'K21495', 'K21494', 'K21493', 'K21492', 'K21491', 'K21490', 'K21489', 'K21488', 'K21487', 'K19687', 'K19168', 'K19165', 'K19164', 'K19163', 'K19161', 'K19159', 'K19158', 'K19156', 'K19155', 'K19094', 'K19093', 'K19092', 'K18923', 'K18918', 'K18843', 'K18842', 'K18840', 'K18839', 'K18831', 'K18830', 'K18829', 'K15773', 'K09159', 'K07746', 'K07723', 'K07334', 'K07172', 'K07062']
			restriction_enzyme = ['K19147', 'K07452', 'K07451', 'K03427', 'K01156', 'K01155', 'K01154', 'K01153']

			for genome in genomes:
				first_viral = 0.0
				check += 1
				if check > 1:
					n -= 14
				all_kegg = 0.0
				all_pfam = 0.0
				all_vog = 0.0
				cat_kegg = 0.0
				cat_pfam = 0.0
				cat_vog = 0.0
				kegg_zero = 0.0
				pfam_zero = 0.0
				pfam_all_counter = 0.0
				pfam_cat_counter = 0.0
				pfam_zero_counter = 0.0
				kegg_all_counter = 0.0
				kegg_cat_counter = 0.0
				kegg_zero_counter = 0.0
				total_genes_counter = 0.0
				annotation_check_counter = 0.0
				kegg_int_rep_count_counter = 0.0
				pfam_int_rep_count_counter = 0.0
				p_v_check_counter = 0.0
				p_k_check_counter = 0.0
				k_v_check_counter = 0.0
				v_check_counter = 0.0
				p_check_counter = 0.0
				k_check_counter = 0.0
				h_check_counter = 0.0
				plasmid_check_counter = 0.0
				redoxin_check_counter = 0.0
				RnR_check_counter = 0.0
				DNA_check_counter = 0.0
				integrase_check_counter = 0.0
				special_check_counter = 0.0
				toxin_check_counter = 0.0
				restriction_check_counter = 0.0
				protein_counter = []
				annotation_check = 0.0
				kegg_int_rep_count = 0.0
				pfam_int_rep_count = 0.0
				p_v_check = 0.0
				p_k_check = 0.0
				k_v_check = 0.0
				v_check = 0.0
				p_check = 0.0
				k_check = 0.0
				h_check = 0.0
				plasmid_check = 0.0
				redoxin_check = 0.0
				RnR_check = 0.0
				DNA_check = 0.0
				integrase_check = 0.0
				special_check = 0.0
				toxin_check = 0.0
				restriction_check = 0.0
				total_genes = 0.0

	####### Number of annotations
				counter = 0
				while counter == 0:
					kegg_check = False
					pfam_check = False
					vog_check = False
					annotation_check_TF = False
					kegg_int_rep_count_TF = False
					pfam_int_rep_count_TF = False
					p_v_check_TF = False
					p_k_check_TF = False
					k_v_check_TF = False
					v_check_TF = False
					p_check_TF = False
					k_check_TF = False
					h_check_TF = False
					plasmid_check_TF = False
					redoxin_check_TF = False
					RnR_check_TF = False
					DNA_check_TF = False
					integrase_check_TF = False
					special_check_TF = False
					toxin_check_TF = False
					restriction_check_TF = False
					if n >= len(annotations):
						break
					if genome == annotations[n]:
						total_genes += 1
						if annotations[n+1] != '':
							kegg_check = True
							all_kegg += 1
							cat_kegg += float(annotations[n+4])
							if annotations[n+4] == "0":
								kegg_zero += 1
							if annotations[n+1] in toxin:
								toxin_check	+= 1
								toxin_check_TF = True
							if annotations[n+1] in restriction_enzyme:
								restriction_check += 1
								restriction_check_TF = True
							if annotations[n+1] in kegg_int_rep_list:
								kegg_int_rep_count += 1
								kegg_int_rep_count_TF = True
						if annotations[n+5] != '':
							pfam_check = True
							all_pfam += 1
							cat_pfam += float(annotations[n+8])
							if annotations[n+8] == "0":
								pfam_zero += 1
							if annotations[n+5] in pfam_int_rep_list:
								pfam_int_rep_count += 1
								pfam_int_rep_count_TF = True
						if annotations[n+9] != '':
							vog_check = True
							all_vog += 1
							cat_vog += float(annotations[n+12])
							if annotations[n+9] in recombinase_transposase:
								plasmid_check += 1
								plasmid_check_TF = True
							if annotations[n+9] in redoxin:
								redoxin_check += 1
								redoxin_check_TF = True
							if annotations[n+9] in RnR:
								RnR_check += 1
								RnR_check_TF = True
							if annotations[n+9] in DNA:
								DNA_check += 1
								DNA_check_TF = True
							if annotations[n+9] in integrase:
								integrase_check += 1
								integrase_check_TF = True
							if annotations[n+9] in special:
								special_check += 1
								special_check_TF = True
						if kegg_check == True and pfam_check == True and vog_check == True:
							annotation_check += 1
							annotation_check_TF = True
						if pfam_check == True and vog_check == True and kegg_check == False:
							p_v_check += 1
							p_v_check_TF = True
						if pfam_check == True and vog_check == False and kegg_check == True:
							p_k_check += 1
							p_k_check_TF = True
						if pfam_check == False and vog_check == True and kegg_check == True:
							k_v_check += 1
							k_v_check_TF = True
						if pfam_check == False and vog_check == True and kegg_check == False:
							v_check += 1
							v_check_TF = True
						if pfam_check == True and vog_check == False and kegg_check == False:
							p_check += 1
							p_check_TF = True
						if pfam_check == False and vog_check == False and kegg_check == True:
							k_check += 1
							k_check_TF = True
						if pfam_check == False and vog_check == False and kegg_check == False:
							h_check += 1
							h_check_TF = True

						if "fragment" in str(genome):
							if first_viral == 0:
								if (kegg_check == True and float(annotations[n+4]) >= 0.1) or (pfam_check == True and float(annotations[n+8]) >= 0.1) or (vog_check == True):
									if len(protein_counter) > 0:
										for item in protein_counter:
											final_proteins.append(item)
									protein_counter = []
									all_kegg += kegg_all_counter
									cat_kegg += kegg_cat_counter
									all_pfam += pfam_all_counter
									cat_pfam += pfam_cat_counter
									kegg_zero += kegg_zero_counter
									pfam_zero += pfam_zero_counter
									total_genes += total_genes_counter
									annotation_check += annotation_check_counter
									kegg_int_rep_count += kegg_int_rep_count_counter
									kegg_int_rep_count += kegg_int_rep_count_counter
									p_v_check += p_v_check_counter
									p_k_check += p_k_check_counter
									k_v_check += k_v_check_counter
									v_check += v_check_counter
									p_check += p_check_counter
									k_check += k_check_counter
									h_check += h_check_counter
									plasmid_check += plasmid_check_counter
									redoxin_check += redoxin_check_counter
									RnR_check += RnR_check_counter
									DNA_check += DNA_check_counter
									integrase_check += integrase_check_counter
									special_check += special_check_counter
									toxin_check += toxin_check_counter
									restriction_check += restriction_check_counter
									pfam_all_counter = 0
									pfam_cat_counter = 0
									kegg_all_counter = 0
									kegg_cat_counter = 0
									pfam_zero_counter = 0
									kegg_zero_counter = 0
									total_genes_counter = 0
									annotation_check_counter = 0.0
									kegg_int_rep_count_counter = 0.0
									pfam_int_rep_count_counter = 0.0
									p_v_check_counter = 0.0
									p_k_check_counter = 0.0
									k_v_check_counter = 0.0
									v_check_counter = 0.0
									p_check_counter = 0.0
									k_check_counter = 0.0
									h_check_counter = 0.0
									plasmid_check_counter = 0.0
									redoxin_check_counter = 0.0
									RnR_check_counter = 0.0
									DNA_check_counter = 0.0
									integrase_check_counter = 0.0
									special_check_counter = 0.0
									total_genes_counter = 0.0
									toxin_check_counter = 0.0
									restriction_check_counter = 0.0
									final_proteins.append(annotations[n-1])
									first_viral = 1
								else:
									all_kegg = 0.0
									all_pfam = 0.0
									cat_kegg = 0.0
									cat_pfam = 0.0
									kegg_zero = 0.0
									pfam_zero = 0.0
									total_genes = 0.0
									annotation_check = 0.0
									kegg_int_rep_count = 0.0
									pfam_int_rep_count = 0.0
									p_v_check = 0.0
									p_k_check = 0.0
									k_v_check = 0.0
									v_check = 0.0
									p_check = 0.0
									k_check = 0.0
									h_check = 0.0
									plasmid_check = 0.0
									redoxin_check = 0.0
									RnR_check = 0.0
									DNA_check = 0.0
									integrase_check = 0.0
									special_check = 0.0
									toxin_check = 0.0
									restriction_check = 0.0

							if first_viral == 1:
								if vog_check == True:
									if len(protein_counter) > 0:
										for item in protein_counter:
											final_proteins.append(item)
									protein_counter = []
									all_kegg += kegg_all_counter
									cat_kegg += kegg_cat_counter
									all_pfam += pfam_all_counter
									cat_pfam += pfam_cat_counter
									kegg_zero += kegg_zero_counter
									pfam_zero += pfam_zero_counter
									total_genes += total_genes_counter
									pfam_all_counter = 0.0
									pfam_cat_counter = 0.0
									kegg_all_counter = 0.0
									kegg_cat_counter = 0.0
									pfam_zero_counter = 0.0
									kegg_zero_counter = 0.0
									total_genes_counter = 0.0
									annotation_check += annotation_check_counter
									kegg_int_rep_count += kegg_int_rep_count_counter
									kegg_int_rep_count += kegg_int_rep_count_counter
									p_v_check += p_v_check_counter
									p_k_check += p_k_check_counter
									k_v_check += k_v_check_counter
									v_check += v_check_counter
									p_check += p_check_counter
									k_check += k_check_counter
									h_check += h_check_counter
									plasmid_check += plasmid_check_counter
									redoxin_check += redoxin_check_counter
									RnR_check += RnR_check_counter
									DNA_check += DNA_check_counter
									integrase_check += integrase_check_counter
									special_check += special_check_counter
									toxin_check += toxin_check_counter
									restriction_check += restriction_check_counter
									pfam_all_counter = 0
									pfam_cat_counter = 0
									kegg_all_counter = 0
									kegg_cat_counter = 0
									pfam_zero_counter = 0
									kegg_zero_counter = 0
									total_genes_counter = 0
									annotation_check_counter = 0.0
									kegg_int_rep_count_counter = 0.0
									pfam_int_rep_count_counter = 0.0
									p_v_check_counter = 0.0
									p_k_check_counter = 0.0
									k_v_check_counter = 0.0
									v_check_counter = 0.0
									p_check_counter = 0.0
									k_check_counter = 0.0
									h_check_counter = 0.0
									plasmid_check_counter = 0.0
									redoxin_check_counter = 0.0
									RnR_check_counter = 0.0
									DNA_check_counter = 0.0
									integrase_check_counter = 0.0
									special_check_counter = 0.0
									total_genes_counter = 0.0
									toxin_check_counter = 0.0
									restriction_check_counter = 0.0
									final_proteins.append(annotations[n-1])
								if vog_check == False:
									if (kegg_check == True and float(annotations[n+4]) >= 0.1) or (pfam_check == True and float(annotations[n+8]) >= 0.1):
										if len(protein_counter) > 0:
											for item in protein_counter:
												final_proteins.append(item)
										protein_counter = []
										all_kegg += kegg_all_counter
										cat_kegg += kegg_cat_counter
										all_pfam += pfam_all_counter
										cat_pfam += pfam_cat_counter
										kegg_zero += kegg_zero_counter
										pfam_zero += pfam_zero_counter
										total_genes += total_genes_counter
										pfam_all_counter = 0.0
										pfam_cat_counter = 0.0
										kegg_all_counter = 0.0
										kegg_cat_counter = 0.0
										pfam_zero_counter = 0.0
										kegg_zero_counter = 0.0
										total_genes_counter = 0.0
										annotation_check += annotation_check_counter
										kegg_int_rep_count += kegg_int_rep_count_counter
										kegg_int_rep_count += kegg_int_rep_count_counter
										p_v_check += p_v_check_counter
										p_k_check += p_k_check_counter
										k_v_check += k_v_check_counter
										v_check += v_check_counter
										p_check += p_check_counter
										k_check += k_check_counter
										h_check += h_check_counter
										plasmid_check += plasmid_check_counter
										redoxin_check += redoxin_check_counter
										RnR_check += RnR_check_counter
										DNA_check += DNA_check_counter
										integrase_check += integrase_check_counter
										special_check += special_check_counter
										toxin_check += toxin_check_counter
										restriction_check += restriction_check_counter
										pfam_all_counter = 0
										pfam_cat_counter = 0
										kegg_all_counter = 0
										kegg_cat_counter = 0
										pfam_zero_counter = 0
										kegg_zero_counter = 0
										total_genes_counter = 0
										annotation_check_counter = 0.0
										kegg_int_rep_count_counter = 0.0
										pfam_int_rep_count_counter = 0.0
										p_v_check_counter = 0.0
										p_k_check_counter = 0.0
										k_v_check_counter = 0.0
										v_check_counter = 0.0
										p_check_counter = 0.0
										k_check_counter = 0.0
										h_check_counter = 0.0
										plasmid_check_counter = 0.0
										redoxin_check_counter = 0.0
										RnR_check_counter = 0.0
										DNA_check_counter = 0.0
										integrase_check_counter = 0.0
										special_check_counter = 0.0
										total_genes_counter = 0.0
										toxin_check_counter = 0.0
										restriction_check_counter = 0.0
										final_proteins.append(annotations[n-1])
									else:
										if kegg_check == True:
											kegg_all_counter += 1
											kegg_cat_counter += float(annotations[n+4])
											all_kegg -= 1
											cat_kegg -= float(annotations[n+4])
											if float(annotations[n+4]) == "0.0":
												kegg_zero -= 1
										if pfam_check == True:
											pfam_all_counter += 1
											pfam_cat_counter += float(annotations[n+8])
											all_pfam -= 1
											cat_pfam -= float(annotations[n+8])
											if float(annotations[n+8]) == "0.0":
												pfam_zero -= 1
										if annotation_check_TF == True:
											annotation_check_counter += 1
											annotation_check -= 1
										if kegg_int_rep_count_TF == True:
											kegg_int_rep_count_counter += 1
											kegg_int_rep_count -= 1
										if pfam_int_rep_count_TF == True:
											pfam_int_rep_count_counter += 1
											pfam_int_rep_count -= 1
										if p_v_check_TF == True:
											p_v_check_counter += 1
											p_v_check -= 1
										if p_k_check_TF == True:
											p_k_check_counter += 1
											p_k_check -= 1
										if k_v_check_TF == True:
											k_v_check_counter += 1
											k_v_check -= 1
										if v_check_TF == True:
											v_check_counter += 1
											v_check -= 1
										if p_check_TF == True:
											p_check_counter += 1
											p_check -= 1
										if k_check_TF == True:
											k_check_counter += 1
											k_check -= 1
										if h_check_TF == True:
											h_check_counter += 1
											h_check -= 1
										if plasmid_check_TF == True:
											plasmid_check_counter += 1
											plasmid_check -= 1
										if redoxin_check_TF == True:
											redoxin_check_counter += 1
											redoxin_check -= 1
										if RnR_check_TF == True:
											RnR_check_counter += 1
											RnR_check -= 1
										if DNA_check_TF == True:
											DNA_check_counter += 1
											DNA_check -= 1
										if integrase_check_TF == True:
											integrase_check_counter += 1
											integrase_check -= 1
										if special_check_TF == True:
											special_check_counter += 1
											special_check -= 1
										if toxin_check_TF == True:
											toxin_check_counter += 1
											toxin_check -= 1
										if restriction_check_TF == True:
											restriction_check_counter += 1
											restriction_check -= 1
										total_genes_counter += 1
										total_genes -= 1
										protein_counter.append(annotations[n-1])

						if "fragment" not in str(genome):
							final_proteins.append(annotations[n-1])
					if genome != annotations[n]:
						counter += 1
					n += 14

	############# Remove contigs that are definitely viral
				if total_genes >= int(orf_low):
					unmod_results.write(str(genome) + '\t' + str(total_genes) + '\t' + str(all_kegg) + '\t' + str(cat_kegg) + '\t' + str(all_pfam) + '\t' + str(cat_pfam) + '\t' + str(all_vog) + '\t' + str(cat_vog) + '\t' + str(kegg_int_rep_count) + '\t' + str(kegg_zero) + '\t' + str(pfam_int_rep_count) + '\t' + str(pfam_zero) + '\t' + str(redoxin_check) + '\t' + str(plasmid_check) + '\t' + str(integrase_check) + '\t' + str(RnR_check) + '\t' + str(DNA_check) + '\t' + str(restriction_check) + '\t' + str(toxin_check) + '\t' + str(special_check) + '\t' + str(annotation_check) +'\t' + str(p_v_check) + '\t' + str(p_k_check) + '\t' + str(k_v_check) + '\t' + str(k_check) + '\t' + str(p_check) + '\t' + str(v_check) + '\t' + str(h_check) + '\n')
	########################### Begin to normalize data #######################
					#else:
					if cat_vog > 0:
						cat_vog = float(cat_vog/all_vog)/float(total_genes)
					if cat_pfam > 0:
						cat_pfam = float(cat_pfam/all_pfam)/float(total_genes)
					if cat_kegg > 0:
						cat_kegg = float(cat_kegg/all_kegg)/float(total_genes)
					all_kegg /= float(total_genes)
					all_pfam /= float(total_genes)
					all_vog /= float(total_genes)
					kegg_zero /= float(total_genes)
					pfam_zero /= float(total_genes)
					annotation_check /= float(total_genes)
					kegg_int_rep_count /= float(total_genes)
					pfam_int_rep_count /= float(total_genes)
					p_v_check /= float(total_genes)
					p_k_check /= float(total_genes)
					k_v_check /= float(total_genes)
					v_check /= float(total_genes)
					p_check /= float(total_genes)
					k_check /= float(total_genes)
					h_check /= float(total_genes)
					plasmid_check /= float(total_genes)
					redoxin_check /= float(total_genes)
					RnR_check /= float(total_genes)
					DNA_check /= float(total_genes)
					integrase_check /= float(total_genes)
					special_check /= float(total_genes)
					toxin_check /= float(total_genes)
					restriction_check /= float(total_genes)
					machine_check += 1
					results.write(str(genome) + '\t' + str(math.log10(total_genes)) + '\t' + str(all_kegg) + '\t' + str(cat_kegg) + '\t' + str(all_pfam) + '\t' + str(cat_pfam) + '\t' + str(all_vog) + '\t' + str(cat_vog) + '\t' + str(kegg_int_rep_count) + '\t' + str(kegg_zero) + '\t' + str(pfam_int_rep_count) + '\t' + str(pfam_zero) + '\t' + str(redoxin_check) + '\t' + str(plasmid_check) + '\t' + str(integrase_check) + '\t' + str(RnR_check) + '\t' + str(DNA_check) + '\t' + str(restriction_check) + '\t' + str(toxin_check) + '\t' + str(special_check) + '\t' + str(annotation_check) +'\t' + str(p_v_check) + '\t' + str(p_k_check) + '\t' + str(k_v_check) + '\t' + str(k_check) + '\t' + str(p_check) + '\t' + str(v_check) + '\t' + str(h_check) + '\n')
					ORFs.append(total_genes)

if len(ORFs) == 0:
	exit()

### Run neural network (MLPClassifier)
plasmid_pred = []
organism_pred = []
if machine_check > 0:
	with open(str(path)+'temp_VIBRANT_results.' + str(base) + '.txt', 'r') as machine_results:
		with open(model, 'rb') as read_model:
			with open(str(path)+'temp_VIBRANT_machine.' + str(base) + '.txt', 'w') as write_file:
				annotation_result = pd.read_csv(machine_results, sep='\t', header=0)
				if len(list(annotations)) > 0:
					anno = preprocessing.normalize(annotation_result[['total genes' , 'all KEGG' , 'category KEGG' , 'all Pfam' , 'category vPfam' , 'all VOG' , 'category VOG' , 'KEGG int-rep' , 'KEGG zero' , 'Pfam int-rep' , 'Pfam zero' , 'VOG redoxin' , 'VOG rec-tran' , 'VOG int' , 'VOG RnR' , 'VOG DNA' , 'KEGG restriction check', 'KEGG toxin check', 'VOG special', 'annotation check' , 'p_v check' , 'p_k check' , 'k_v check' , 'k check' , 'p check' , 'v check' , 'h check']].values)
					anno = pd.DataFrame(anno, columns=['total genes' , 'all KEGG' , 'category KEGG' , 'all Pfam' , 'category vPfam' , 'all VOG' , 'category VOG' , 'KEGG int-rep' , 'KEGG zero' , 'Pfam int-rep' , 'Pfam zero' , 'VOG redoxin' , 'VOG rec-tran' , 'VOG int' , 'VOG RnR' , 'VOG DNA' , 'KEGG restriction check', 'KEGG toxin check', 'VOG special', 'annotation check' , 'p_v check' , 'p_k check' , 'k_v check' , 'k check' , 'p check' , 'v check' , 'h check'])
					anno[['genome']] = annotation_result[['genome']]
					anno_run = anno.drop(['genome'], axis='columns')
					anno_pred = anno['genome']
					model = pickle.load(read_model)
					predictions = model.predict(anno_run)
				else:
					predictions = []
				n = 0
				if len(predictions) > 0:
					for item in predictions:
						write_file.write(str(anno_pred[n]).replace("$~&", " ").replace('^@%','"') + "\t" + str(predictions[n]) + "\n")
						n += 1
				else:
					write_file.write('')

	with open(str(path)+'temp_VIBRANT_machine.' + str(base) + '.txt', 'r') as machine_file:
		machine = machine_file.read().replace('\n','\t').split('\t')
		n = 1
		while n < len(machine):
			if machine[n].replace(" ", "$~&").replace('"','^@%') == 'virus':
				virus_pred.append(str(machine[n-1]).replace(" ", "$~&").replace('"','^@%'))
			elif machine[n].replace(" ", "$~&").replace('"','^@%') == 'plasmid':
				plasmid_pred.append(str(machine[n-1]).replace(" ", "$~&").replace('"','^@%'))
			elif machine[n].replace(" ", "$~&").replace('"','^@%') == 'organism':
				organism_pred.append(str(machine[n-1]).replace(" ", "$~&").replace('"','^@%'))
			n += 2


final_check = []
for item in temp_list:
	if str(item) in final_proteins:
		final_check.append(item)


with open(str(path)+'temp1_VIBRANT_annotations.' + str(base) + '.txt', 'r') as temp:
	with open(str(path)+'temp2_VIBRANT_annotations.' + str(base) + '.txt', 'w') as final_anno:
		for item in temp_list:
			if str(item) in final_check:
				final_anno.write(str(item) + "\t" + str(temp_dict[item]) + "\n")

with open(str(path)+'temp2_VIBRANT_annotations.' + str(base) + '.txt', 'r') as final_anno:
	annotations = final_anno.read().replace('\n','\t').split('\t')
	annotations = annotations[:-1]

subprocess.run('rm '+str(path)+'temp1_VIBRANT_annotations.' + str(base) + '.txt', shell=True)

n = 1
i = 0
full_genomes = list()
genomes = list()
while n < len(annotations):
	full_genomes.append(annotations[n])
	n += 14
for item in full_genomes:
	if item not in genomes:
		genomes.append(item)
n = 1
check = 0
viral_genomes = []
viral_proteins = []
final = []

integrase = ['VOG15133', 'VOG20969', 'VOG02658', 'VOG04024', 'VOG01778', 'VOG02371','VOG00041']
kegg_int_rep_list = ['K07729', 'K03733', 'K04763', 'K14059', 'K21039', 'K01356', 'K18918', 'K07741', 'K21528', 'K06400', 'K01152', 'K07481', 'K07482', 'K07483', 'K07484', 'K07485', 'K07486', 'K07487', 'K07488', 'K07489', 'K07491', 'K07492', 'K07493', 'K07494', 'K07495', 'K07496', 'K07497', 'K07498', 'K07499', 'K18320', 'K23209']
pfam_int_rep_list = ['PF00239.21', 'PF00589.22', 'PF00665.26', 'PF00872.18', 'PF01385.19', 'PF01526.17', 'PF01527.20', 'PF01548.17', 'PF01609.21', 'PF01610.17', 'PF01726.16', 'PF01797.16', 'PF01818.17', 'PF02061.16', 'PF02316.16', 'PF02899.17', 'PF02914.15', 'PF03050.14', 'PF03374.14', 'PF03400.13', 'PF04754.12', 'PF04761.12', 'PF05269.11', 'PF05598.11', 'PF05973.14', 'PF06543.12', 'PF07022.13', 'PF07282.11', 'PF07508.13', 'PF09003.10', 'PF09299.11', 'PF09588.10', 'PF09669.10', 'PF10551.9', 'PF12167.8', 'PF12472.8', 'PF12760.7', 'PF13022.6', 'PF13333.6', 'PF13340.6', 'PF13356.6', 'PF13408.6', 'PF13586.6', 'PF13701.6', 'PF13737.6', 'PF14657.6', 'PF14659.6', 'PF16452.5', 'PF16795.5']
DNA = ['VOG02853', 'VOG07917', 'VOG24666', 'VOG03664', 'VOG04691', 'VOG21077', 'VOG00821', 'VOG00163', 'VOG11329', 'VOG11335', 'VOG10269', 'VOG11468', 'VOG00817', 'VOG01029', 'VOG00031', 'VOG01431', 'VOG03850', 'VOG00038', 'VOG00108', 'VOG01695', 'VOG00654', 'VOG01045', 'VOG00404', 'VOG00029', 'VOG00056', 'VOG00079', 'VOG00084', 'VOG00109', 'VOG00272', 'VOG00329', 'VOG00356', 'VOG00422', 'VOG00474', 'VOG00476', 'VOG00532', 'VOG00614', 'VOG00674', 'VOG00709', 'VOG00724', 'VOG00728', 'VOG00815', 'VOG00863', 'VOG00895', 'VOG00908', 'VOG00935', 'VOG01042', 'VOG01137', 'VOG01374', 'VOG01384', 'VOG01422', 'VOG01450', 'VOG01519', 'VOG01585', 'VOG01605', 'VOG01637', 'VOG01664', 'VOG01697', 'VOG01763', 'VOG01852', 'VOG01895', 'VOG01901', 'VOG02015', 'VOG02108', 'VOG02209', 'VOG02395', 'VOG02401', 'VOG02569', 'VOG02821', 'VOG02851', 'VOG02874', 'VOG03009', 'VOG03044', 'VOG03307', 'VOG03317', 'VOG03439', 'VOG03479', 'VOG03586', 'VOG03953', 'VOG04002', 'VOG04063', 'VOG04302', 'VOG04652', 'VOG04729', 'VOG04850', 'VOG04962', 'VOG04998', 'VOG05011', 'VOG05034', 'VOG05091', 'VOG05095', 'VOG05100', 'VOG05105', 'VOG05115', 'VOG05119', 'VOG05300', 'VOG05317', 'VOG05323', 'VOG05344', 'VOG05382', 'VOG05406', 'VOG05693', 'VOG05802', 'VOG05966', 'VOG05972', 'VOG06283', 'VOG07014', 'VOG07060', 'VOG07770', 'VOG08033', 'VOG08045', 'VOG08144', 'VOG08273', 'VOG08380', 'VOG08409', 'VOG08440', 'VOG08827', 'VOG09260', 'VOG09702', 'VOG09774', 'VOG10038', 'VOG10089', 'VOG10101', 'VOG10140', 'VOG10207', 'VOG10235', 'VOG10289', 'VOG11125', 'VOG11285', 'VOG11351', 'VOG11729', 'VOG11958', 'VOG12206', 'VOG12269', 'VOG14398', 'VOG14781', 'VOG17046', 'VOG17463', 'VOG17529', 'VOG18148', 'VOG18226', 'VOG18480', 'VOG18481', 'VOG18702', 'VOG18819', 'VOG18988', 'VOG19293', 'VOG19373', 'VOG19571', 'VOG19864', 'VOG20018', 'VOG20102', 'VOG20254', 'VOG20418', 'VOG20985', 'VOG21265', 'VOG21353', 'VOG21722', 'VOG21964', 'VOG22116', 'VOG22861', 'VOG22863', 'VOG23467', 'VOG23468', 'VOG23717', 'VOG24749', 'VOG24750', 'VOG24981', 'VOG00393', 'VOG24404', 'VOG10258', 'VOG00557', 'VOG02455', 'VOG03269', 'VOG06969', 'VOG03619', 'VOG12136', 'VOG19468', 'VOG04560', 'VOG15690', 'VOG00222', 'VOG00157', 'VOG01469', 'VOG02646', 'VOG12013', 'VOG08321', 'VOG12635', 'VOG03573', 'VOG00140', 'VOG16387', 'VOG11378', 'VOG24766', 'VOG10041', 'VOG09778', 'VOG03513', 'VOG09828', 'VOG01256', 'VOG01464', 'VOG02533', 'VOG05391', 'VOG07031', 'VOG00098']
special = ['VOG17265', 'VOG23378', 'VOG03971', 'VOG19974', 'VOG21494', 'VOG19846', 'VOG21010', 'VOG23255', 'VOG23256', 'VOG08837', 'VOG08834', 'VOG21336', 'VOG19582', 'VOG16059', 'VOG16060', 'VOG03412', 'VOG02335', 'VOG23390', 'VOG24774', 'VOG24789', 'VOG23939', 'VOG19698', 'VOG20700', 'VOG21766', 'VOG24620', 'VOG21782', 'VOG24791', 'VOG23557', 'VOG11305', 'VOG22164', 'VOG02628', 'VOG15962', 'VOG17611', 'VOG18579', 'VOG22447', 'VOG23134', 'VOG10066', 'VOG04780', 'VOG18550', 'VOG18551', 'VOG23079', 'VOG23391', 'VOG24118', 'VOG21998', 'VOG06315', 'VOG10278', 'VOG16253', 'VOG19887', 'VOG21123', 'VOG02491', 'VOG19122', 'VOG03884', 'VOG15965', 'VOG11324', 'VOG22086', 'VOG04393', 'VOG19193', 'VOG20183', 'VOG04197', 'VOG19981', 'VOG14348', 'VOG22444', 'VOG04150', 'VOG18810', 'VOG11983', 'VOG12138', 'VOG10279', 'VOG19484', 'VOG23315', 'VOG23490', 'VOG00528', 'VOG16358', 'VOG04297', 'VOG04356', 'VOG23566', 'VOG22745', 'VOG03721', 'VOG11306', 'VOG19832', 'VOG17540', 'VOG00955', 'VOG04139', 'VOG06441', 'VOG06443', 'VOG07385', 'VOG07388', 'VOG07473', 'VOG08323', 'VOG08324', 'VOG08325', 'VOG08326', 'VOG08327', 'VOG08328', 'VOG08340', 'VOG08343', 'VOG14347', 'VOG15036', 'VOG15928', 'VOG15935', 'VOG15939', 'VOG15963', 'VOG19975', 'VOG21139', 'VOG22650', 'VOG12902', 'VOG01051', 'VOG19855', 'VOG11420', 'VOG15930', 'VOG08329', 'VOG02389', 'VOG01861', 'VOG01052', 'VOG04138', 'VOG15936', 'VOG08336', 'VOG15934', 'VOG08865', 'VOG21079', 'VOG23556', 'VOG08341', 'VOG21085', 'VOG15964', 'VOG23077', 'VOG15933', 'VOG15938', 'VOG08332', 'VOG08335', 'VOG02193', 'VOG16504', 'VOG21081', 'VOG04415', 'VOG09120', 'VOG21080', 'VOG08333', 'VOG21084', 'VOG21086', 'VOG04907', 'VOG08334', 'VOG04822', 'VOG04228', 'VOG21088', 'VOG03844', 'VOG04764', 'VOG21082', 'VOG04935', 'VOG03204', 'VOG01734', 'VOG00296', 'VOG03147', 'VOG06128', 'VOG06289', 'VOG17834', 'VOG23226', 'VOG13979', 'VOG18491', 'VOG24613', 'VOG06316', 'VOG06375', 'VOG06213', 'VOG07922', 'VOG02468', 'VOG21446', 'VOG12267', 'VOG03194', 'VOG01505', 'VOG06374', 'VOG06140', 'VOG15391', 'VOG03520', 'VOG18860', 'VOG19382', 'VOG11975', 'VOG06137', 'VOG15761', 'VOG02145', 'VOG04079', 'VOG07639', 'VOG06285', 'VOG12899', 'VOG23376', 'VOG24353', 'VOG18667', 'VOG00367', 'VOG00227', 'VOG19851', 'VOG23555', 'VOG12225', 'VOG23374', 'VOG01647', 'VOG22306', 'VOG21440', 'VOG22879', 'VOG20648', 'VOG24612', 'VOG01694', 'VOG11430', 'VOG11302', 'VOG03777', 'VOG13935', 'VOG04522', 'VOG10851', 'VOG00992', 'VOG17613', 'VOG02117', 'VOG03059', 'VOG03404', 'VOG03926', 'VOG19265', 'VOG19587', 'VOG23227', 'VOG24068', 'VOG23257', 'VOG01346', 'VOG03399', 'VOG01179', 'VOG01973', 'VOG19429', 'VOG05705', 'VOG09906', 'VOG11707', 'VOG21325', 'VOG10212', 'VOG06218', 'VOG11404', 'VOG00156', 'VOG00293', 'VOG00541', 'VOG00564', 'VOG01022', 'VOG01143', 'VOG01447', 'VOG01481', 'VOG01933', 'VOG01941', 'VOG02050', 'VOG02069', 'VOG02360', 'VOG02603', 'VOG02784', 'VOG03114', 'VOG03170', 'VOG03254', 'VOG03303', 'VOG03534', 'VOG04291', 'VOG12621', 'VOG14566', 'VOG14730', 'VOG16449', 'VOG17449', 'VOG17452', 'VOG18314', 'VOG19081', 'VOG19271', 'VOG19806', 'VOG20055', 'VOG20145', 'VOG20156', 'VOG21749', 'VOG21991', 'VOG22176', 'VOG22825', 'VOG23189', 'VOG23575', 'VOG23576', 'VOG23927', 'VOG24826', 'VOG02754', 'VOG02419', 'VOG01341', 'VOG03084', 'VOG24311', 'VOG11764', 'VOG18132', 'VOG23129', 'VOG16767', 'VOG01035', 'VOG17546', 'VOG00257', 'VOG19860', 'VOG02944', 'VOG17022', 'VOG18956', 'VOG00101', 'VOG04298', 'VOG00651', 'VOG17451', 'VOG12616', 'VOG04398', 'VOG12622', 'VOG01720', 'VOG03069', 'VOG00929', 'VOG01372', 'VOG11405', 'VOG03383', 'VOG18346', 'VOG03470', 'VOG23090', 'VOG11456', 'VOG11458', 'VOG03487', 'VOG03143', 'VOG11457', 'VOG22335', 'VOG18264', 'VOG04815', 'VOG12907', 'VOG12072', 'VOG01537', 'VOG13937', 'VOG01910', 'VOG01795', 'VOG02823', 'VOG11683', 'VOG21902', 'VOG02067', 'VOG19992', 'VOG01241', 'VOG12227', 'VOG09900', 'VOG12228', 'VOG11088', 'VOG11982', 'VOG04898', 'VOG19579', 'VOG09336', 'VOG01975', 'VOG07839', 'VOG15568', 'VOG13642', 'VOG15567', 'VOG21165', 'VOG02813', 'VOG03014', 'VOG03217', 'VOG10214', 'VOG10242', 'VOG17242', 'VOG23908', 'VOG10146', 'VOG00489', 'VOG01633', 'VOG09941', 'VOG11985', 'VOG17915', 'VOG02782', 'VOG21767', 'VOG17917', 'VOG01513', 'VOG21669', 'VOG10251', 'VOG10252', 'VOG12917', 'VOG10172', 'VOG17900', 'VOG06155', 'VOG01389', 'VOG10243', 'VOG11711', 'VOG17926', 'VOG15597', 'VOG17914', 'VOG12918', 'VOG23907', 'VOG02413', 'VOG01350', 'VOG02080', 'VOG01771', 'VOG10174', 'VOG00501', 'VOG18043', 'VOG23031', 'VOG23235', 'VOG23375', 'VOG13980', 'VOG18538', 'VOG21509', 'VOG17756', 'VOG10185', 'VOG13943', 'VOG15967', 'VOG18134', 'VOG10265', 'VOG14735', 'VOG24473', 'VOG15541', 'VOG03699', 'VOG24632', 'VOG14594', 'VOG13649', 'VOG10325', 'VOG10538', 'VOG12146', 'VOG19196', 'VOG19985', 'VOG10280', 'VOG22425', 'VOG19483', 'VOG18036', 'VOG12921', 'VOG16362', 'VOG10264', 'VOG01434', 'VOG13459', 'VOG22145', 'VOG21505', 'VOG02996', 'VOG16532', 'VOG04492', 'VOG20016', 'VOG20963', 'VOG22173', 'VOG24447', 'VOG00549', 'VOG11341', 'VOG11344', 'VOG15577', 'VOG23471', 'VOG23462', 'VOG12900', 'VOG20627', 'VOG10260', 'VOG10188', 'VOG00765', 'VOG10147', 'VOG12432', 'VOG01376', 'VOG01809', 'VOG00090', 'VOG03224', 'VOG02684', 'VOG01873', 'VOG01543', 'VOG23273', 'VOG12190', 'VOG10104', 'VOG03310', 'VOG03354', 'VOG01249', 'VOG11459', 'VOG01806', 'VOG11223', 'VOG10106', 'VOG11579', 'VOG11581', 'VOG02473', 'VOG01629', 'VOG02333', 'VOG11583', 'VOG11585', 'VOG02914', 'VOG08331', 'VOG03803', 'VOG10263', 'VOG10262', 'VOG01105', 'VOG11582', 'VOG10261', 'VOG15579', 'VOG02539', 'VOG01321', 'VOG19840', 'VOG15362', 'VOG04533', 'VOG02548', 'VOG20766', 'VOG03436', 'VOG01591', 'VOG11765', 'VOG00959', 'VOG11303', 'VOG24545', 'VOG24825', 'VOG23935', 'VOG02701', 'VOG02490', 'VOG01437', 'VOG01551', 'VOG11142', 'VOG18035', 'VOG00319', 'VOG01793', 'VOG11125', 'VOG11089', 'VOG00067', 'VOG10266', 'VOG00637', 'VOG00748', 'VOG00988', 'VOG00639', 'VOG00364', 'VOG01829', 'VOG01954', 'VOG00095', 'VOG01944', 'VOG01489', 'VOG12605', 'VOG00678', 'VOG04497', 'VOG00192', 'VOG13487', 'VOG01567', 'VOG02076', 'VOG01483', 'VOG03356', 'VOG00700', 'VOG00720', 'VOG00922', 'VOG00167', 'VOG14056', 'VOG14059', 'VOG14055', 'VOG14054', 'VOG01405', 'VOG14051', 'VOG03000', 'VOG14061', 'VOG00710', 'VOG00547', 'VOG00359', 'VOG02998', 'VOG01571', 'VOG02790', 'VOG23481', 'VOG01003', 'VOG11749', 'VOG19498', 'VOG02752', 'VOG00871', 'VOG01495', 'VOG23339', 'VOG00046', 'VOG00530', 'VOG00177', 'VOG10015', 'VOG09766', 'VOG11785', 'VOG01021', 'VOG01027', 'VOG10954', 'VOG10953', 'VOG00460', 'VOG04756', 'VOG04499', 'VOG00681', 'VOG01398', 'VOG10478', 'VOG10026', 'VOG01666', 'VOG02486', 'VOG00106', 'VOG02984', 'VOG00923', 'VOG02623', 'VOG00392', 'VOG00534', 'VOG02574', 'VOG09790', 'VOG19369', 'VOG09812', 'VOG10752', 'VOG02457', 'VOG09829', 'VOG10761', 'VOG00132', 'VOG00334', 'VOG00186', 'VOG00360', 'VOG02510', 'VOG01904', 'VOG09805', 'VOG09806', 'VOG01488', 'VOG02498', 'VOG00721', 'VOG03347', 'VOG01104', 'VOG01630', 'VOG03107', 'VOG10125', 'VOG01889', 'VOG10754', 'VOG10758', 'VOG02435', 'VOG00419', 'VOG00608', 'VOG10316', 'VOG03023', 'VOG00644', 'VOG10950', 'VOG01039', 'VOG03428', 'VOG00179', 'VOG00305', 'VOG00743', 'VOG01162', 'VOG02913', 'VOG03429', 'VOG01797', 'VOG01359', 'VOG11082', 'VOG00218', 'VOG00149', 'VOG10315', 'VOG03272', 'VOG00951', 'VOG00699', 'VOG11949', 'VOG00283', 'VOG00698', 'VOG11950', 'VOG00410', 'VOG00793', 'VOG00097', 'VOG02938', 'VOG01397', 'VOG01966', 'VOG03252', 'VOG10951', 'VOG00230', 'VOG01803', 'VOG02869', 'VOG02762', 'VOG00912', 'VOG00129', 'VOG10317', 'VOG03064', 'VOG00125', 'VOG00381', 'VOG00687', 'VOG00751', 'VOG03528', 'VOG02524', 'VOG09814', 'VOG00744', 'VOG01880', 'VOG00537', 'VOG00686', 'VOG00991', 'VOG09809', 'VOG00066', 'VOG01144', 'VOG00819', 'VOG10319', 'VOG00158', 'VOG02516', 'VOG00771', 'VOG00295', 'VOG00981', 'VOG02615', 'VOG00069', 'VOG02877', 'VOG00243', 'VOG03165', 'VOG00324', 'VOG10433', 'VOG00241', 'VOG00308', 'VOG00677', 'VOG00326', 'VOG00727', 'VOG00470', 'VOG01125', 'VOG00236', 'VOG00633', 'VOG01646', 'VOG00213', 'VOG00064', 'VOG09810', 'VOG01297', 'VOG00047', 'VOG00921', 'VOG00321', 'VOG01965', 'VOG01802', 'VOG04946']
toxin = ['K21498', 'K21497', 'K21496', 'K21495', 'K21494', 'K21493', 'K21492', 'K21491', 'K21490', 'K21489', 'K21488', 'K21487', 'K19687', 'K19168', 'K19165', 'K19164', 'K19163', 'K19161', 'K19159', 'K19158', 'K19156', 'K19155', 'K19094', 'K19093', 'K19092', 'K18923', 'K18918', 'K18843', 'K18842', 'K18840', 'K18839', 'K18831', 'K18830', 'K18829', 'K15773', 'K09159', 'K07746', 'K07723', 'K07334', 'K07172', 'K07062']
restriction_enzyme = ['K19147', 'K07452', 'K07451', 'K03427', 'K01156', 'K01155', 'K01154', 'K01153']

prophage_integrase = []

with open(infile + '_genome_quality.out', 'w') as quality:
	with open(infile + '.phages_combined.out', 'w') as accnos:
		for genome in genomes:
			check += 1
			if check > 1:
				n -= 14
			kegg_zero = 0
			pfam_zero = 0
			all_kegg = 0
			all_pfam = 0
			all_vog = 0
			cat_kegg = 0
			cat_vog = 0
			cat_pfam = 0
			counter = 0
			DNA_check = 0
			total_genes = 0
			annotation_check = 0
			v_check = 0
			h_check = 0
			k_v_check = 0
			kegg_int_rep_count = 0
			pfam_int_rep_count = 0
			vog_int_count = 0
			toxin_check = 0
			restriction_check = 0
			special_check = 0
			confidence = False

			while counter == 0:
				kegg_check = False
				pfam_check = False
				vog_check = False
				if n >= len(annotations):
					break
				if genome == annotations[n]:
					total_genes += 1
					if annotations[n+1] != '':
						all_kegg += 1
						kegg_check = True
						if annotations[n+4] == "0":
							kegg_zero += 1
						cat_kegg += float(annotations[n+4])
						if annotations[n+1] in kegg_int_rep_list:
							kegg_int_rep_count += 1
						if annotations[n+1] in toxin:
							toxin_check += 1
						if annotations[n+1] in restriction_enzyme:
							restriction_check += 1
					if annotations[n+5] != '':
						all_pfam += 1
						pfam_check = True
						if annotations[n+8] == "0":
							pfam_zero += 1
						cat_pfam += float(annotations[n+8])
						if annotations[n+5] in pfam_int_rep_list:
							pfam_int_rep_count += 1
					if annotations[n+9] != '':
						all_vog += 1
						vog_check = True
						cat_vog += float(annotations[n+12])
						if annotations[n+9] in integrase:
							vog_int_count += 1
							if genome not in prophage_integrase:
								prophage_integrase.append(genome)
						if annotations[n+9] in DNA:
							DNA_check += 1
						if annotations[n+9] in special:
							special_check += 1
					if kegg_check == True and pfam_check == True and vog_check == True:
						annotation_check += 1
					if vog_check == True and kegg_check == False and pfam_check == False:
						v_check += 1
					if vog_check == False and kegg_check == False and pfam_check == False:
						h_check += 1
					if vog_check == True and kegg_check == True and pfam_check == False:
						k_v_check += 1
				if genome != annotations[n]:
					counter += 1

				n += 14

			if virome == False:
				if genome in virus_pred:
					confidence = True
					if all_vog > 3:
						if all_kegg > all_vog and all_pfam > all_vog and annotation_check/all_vog >= 0.9 and cat_kegg < 0.1 and cat_pfam < 0.1:
							confidence = False
						if special_check == 0 and (toxin_check >= 2 or restriction_check >= 2 or toxin_check+restriction_check >= 3) and annotation_check/all_vog >= 0.9 and all_kegg > all_vog and all_pfam > all_vog:
							confidence = False
						if kegg_zero > all_vog and pfam_zero > all_vog and total_genes >= 10:
							confidence = False
						if annotation_check/all_vog >= 0.9 and all_kegg > all_vog and all_pfam > all_vog and h_check/total_genes < 0.2 and total_genes >= 10:
							confidence = False
					if all_vog == 0 and all_kegg == 0 and all_pfam == 0:
						confidence = False
					if (all_vog == 0 or (cat_vog == 0 and all_vog <= 1)) and (all_kegg == 0 or all_kegg < 3) and (all_pfam == 0 or all_pfam < 3) and all_kegg > all_vog and all_pfam > all_vog:
						confidence = False
					if all_vog == 0 and total_genes == 4:
						confidence = False
					if (pfam_int_rep_count > 0 or kegg_int_rep_count > 0) and vog_int_count == 0 and all_vog < all_kegg and all_vog < all_pfam and special_check == 0 and cat_kegg < 2 and cat_pfam < 2:
						confidence = False
					if total_genes <= 5 and all_vog <= all_kegg or all_vog <= all_pfam and DNA_check == all_vog:
						confidence = False

				if genome in plasmid_pred:
					if total_genes >= 20 and vog_int_count >= 1 and all_vog > all_kegg and (all_vog > all_pfam or cat_pfam > 15) and cat_vog > 20 and (v_check >= 3 or cat_vog > 80):
						confidence = True
					if vog_int_count >= 1 and all_vog > all_kegg and all_vog > all_pfam and all_vog >= 10:
						confidence = True

				if genome in organism_pred:
					if total_genes >= 30 and all_vog >= 10 and all_vog > all_kegg*2 and all_vog > all_pfam and cat_pfam > 20 and special_check >= 5:
						confidence = True
					if total_genes >= 30 and all_vog >= 10 and h_check/total_genes >= 0.20 and cat_pfam > 10 and cat_vog > 30 and cat_kegg > 8 and v_check > 1:
						confidence = True
					if total_genes >= 30 and all_vog >= 10 and all_vog > all_kegg and cat_kegg >= 10 and cat_pfam >= 10 and v_check > 5:
						confidence = True
					if total_genes >= 10 and all_vog >= 8 and all_vog >= all_kegg*1.5 and all_vog >= all_pfam*1.5 and cat_kegg >= 1 and cat_pfam >= 1:
						confidence = True
					if total_genes >= 30 and v_check/total_genes >= 0.25 and all_vog > all_kegg and all_vog > all_pfam:
						confidence = True
					if total_genes >= 10 and special_check > 0 and DNA_check > 0 and all_vog > all_kegg and all_vog > all_pfam and v_check > 0 and annotation_check < all_vog:
						confidence = True
					if h_check/total_genes >= 0.75 and total_genes > 20 and all_vog > all_kegg and all_vog > all_pfam and cat_kegg > 1 and cat_pfam > 1:
						confidence = True
					if all_vog == 0 and total_genes == 4 and cat_pfam < 20 and cat_kegg < 20:
						confidence = False

				if (kegg_zero > 10 or pfam_zero > 1) and (all_pfam >= all_vog*1.5 or all_kegg >= all_vog*1.5):
					confidence = False
				if all_vog > 10 and all_pfam >= all_vog*2 or all_kegg >= all_vog*2:
					confidence = False
				if all_vog > 10 and all_pfam > all_vog*1.5 and all_kegg > all_vog*1.5:
					confidence = False

			elif virome == True:
				if genome in virus_pred:
					confidence = True
				if genome in plasmid_pred:
					confidence = True
					if kegg_zero >= 10 or pfam_zero >= 10:
						confidence = False
					if (all_kegg > 5 and all_kegg > all_vog*3) or (all_pfam > 5 and all_pfam > all_vog*3):
						confidence = False
				elif genome in organism_pred:
					if all_vog >= all_kegg-1 and all_vog >= all_pfam-1:
						confidence = True
					elif all_vog == 0 and all_kegg == 0 and all_pfam == 0:
						confidence = True
					elif all_vog > 0 and all_kegg/total_genes < 0.1 and all_kegg/total_genes < 0.1 and (h_check/total_genes >= 0.5 or all_vog/total_genes >= 0.2):
						confidence = True
					elif all_vog >= all_kegg-1 and all_vog >= all_pfam-1 and cat_kegg >= 1 and cat_pfam >= 1:
						confidence = True
					elif total_genes >= 30 and v_check/total_genes >= 0.1 and all_vog >= all_kegg-2 and all_vog >= all_pfam-2:
						confidence = True
					elif special_check >= 2 and v_check >= 2:
						confidence = True
					elif h_check/total_genes >= 0.6 and total_genes > 10 and (cat_kegg > 1 or all_vog >= all_kegg-1) and (cat_pfam > 1 or all_vog >= all_pfam-1):
						confidence = True
					elif all_vog/total_genes >= 0.75 and total_genes >= 10:
						confidence = True
			if "fragment" in str(genome) and total_genes < 8:
				confidence = False

			if confidence == False:
				if genome in prophage_integrase:
					prophage_integrase.remove(genome)
			if confidence == True:
				viral_genomes.append(genome)
				accnos.write(str(genome).replace("$~&", " ").replace('^@%','"') + "\n")
				style = 'lytic'
				if genome in prophage_integrase or "fragment" in str(genome):
					style = 'lysogenic'
				if DNA_check >= 5 and special_check >= 10 and total_genes >= 30 and all_vog >= 10 and v_check >= 5:
					quality.write(str(genome.replace("$~&", " ").replace('^@%','"')) + "\t" + str(style) + "\t" + "high quality draft\n")
				elif special_check >= 4 and DNA_check >= 2 and all_vog >= 10 and v_check >= 1:
					quality.write(str(genome.replace("$~&", " ").replace('^@%','"')) + "\t" + str(style) + "\t" + "high quality draft\n")
				elif DNA_check >= 2 and special_check >= 2 and total_genes >= 10 and all_vog >= 5:
					quality.write(str(genome.replace("$~&", " ").replace('^@%','"')) + "\t" + str(style) + "\t" + "medium quality draft\n")
				elif DNA_check >= 1 and special_check >= 1 and total_genes >= 15 and all_vog >= 5 and v_check >= 2:
					quality.write(str(genome.replace("$~&", " ").replace('^@%','"')) + "\t" + str(style) + "\t" + "medium quality draft\n")
				else:
					quality.write(str(genome.replace("$~&", " ").replace('^@%','"')) + "\t" + str(style) + "\t" + "low quality draft\n")


if len(viral_genomes) == 0:
	subprocess.call('rm '+str(path)+'temp2_VIBRANT_annotations.' + str(base) + '.txt'+' 2>/dev/null', shell=True)
	subprocess.call('rm '+str(path)+'temp_VIBRANT_results.' + str(base) + '.txt'+' 2>/dev/null', shell=True)
	subprocess.call('rm '+str(path)+str(base)+'.pass.faa'+' 2>/dev/null', shell=True)
	subprocess.call('rm '+str(path)+str(base)+'.appended*faa'+' 2>/dev/null', shell=True)
	subprocess.call('rm '+str(path)+str(base)+'.master.txt'+' 2>/dev/null', shell=True)
	exit()

for protein in final_check:
	if protein.rsplit("_",1)[0] in viral_genomes:
		final.append(protein)


with open(AMG_list, 'r') as AMG_file:
	AMG_list = AMG_file.read().split('\n')


with open(str(path)+'temp2_VIBRANT_annotations.' + str(base) + '.txt', 'r') as temp:
	with open(str(path)+'VIBRANT_annotations.' + str(base) + '.txt', 'w') as output:
		with open(annotation_names, 'r') as annotation_names:
			with open(str(path)+'VIBRANT_AMGs.' + str(base) + '.txt', 'w') as AMG_out:
				with open(str(path)+'VIBRANT_genbank_table.' + str(base) + '.txt', 'w') as gb_out:
					anno_dict = {}
					annotation_names = annotation_names.read().replace('\n','\t').split('\t')
					annotation_dict = {}
					AMG_dict = {}
					genbank_dict_full = {}
					genbank_dict = {}
					n = 0
					while n < len(annotation_names):
						annotation_dict.update({annotation_names[n]:annotation_names[n+1]})
						n += 2
					n = 0
					while n < len(annotations):
						AMG = ''
						if annotations[n+2] != '':
							ko_name = annotation_dict[annotations[n+2]]
						else:
							ko_name = ''
						if annotations[n+6] != '':
							pfam_name = annotation_dict[annotations[n+6]]
						else:
							pfam_name = ''
						if annotations[n+10] != '':
							vog_name = annotation_dict[annotations[n+10]]
						else:
							vog_name = ''
						if annotations[n+2] != '':
							if annotations[n+2] in AMG_list:
								AMG = 'AMG'
								AMG_dict.update({str(annotations[n]):str(annotations[n+1]).replace("$~&", " ").replace('^@%','"') + '\t' + str(annotations[n+2]) + '\t' + str(ko_name) + '\t' + str(annotations[n+6]) + '\t' + str(pfam_name)})

						if annotations[n+2] != '':
							genbank_dict_full.update({str(annotations[n]):str(annotations[n+1].replace("$~&", " ").replace('^@%','"'))+"\t"+str(annotations[n+2])+"\t"+str(ko_name)})
							genbank_dict.update({str(annotations[n]):str(annotations[n+2])+"\t"+str(ko_name)})
						elif annotations[n+10] != '':
							genbank_dict_full.update({str(annotations[n]):str(annotations[n+1].replace("$~&", " ").replace('^@%','"'))+"\t"+str(annotations[n+10])+"\t"+str(vog_name)})
							genbank_dict.update({str(annotations[n]):str(annotations[n+10])+"\t"+str(vog_name)})
						elif annotations[n+6] != '':
							genbank_dict_full.update({str(annotations[n]):str(annotations[n+1].replace("$~&", " ").replace('^@%','"'))+"\t"+str(annotations[n+6])+"\t"+str(pfam_name)})
							genbank_dict.update({str(annotations[n]):str(annotations[n+6])+"\t"+str(pfam_name)})
						else:
							genbank_dict_full.update({str(annotations[n]):str(annotations[n+1].replace("$~&", " ").replace('^@%','"'))+"\t"+"None"+"\t"+"hypothetical protein"})
							genbank_dict.update({str(annotations[n]):"None"+"\t"+"hypothetical protein"})

						temp_result = str(annotations[n+1]).replace("$~&", " ").replace('^@%','"') + "\t" + str(annotations[n+2]) + "\t" + str(AMG) + "\t" + str(ko_name) + "\t" + str(annotations[n+3]) + "\t" + str(annotations[n+4]) + "\t" + str(annotations[n+5]) + "\t" + str(annotations[n+6]) + "\t" + str(pfam_name) + "\t" + str(annotations[n+7]) + "\t" + str(annotations[n+8]) + "\t" + str(annotations[n+9]) + "\t" + str(annotations[n+10]) + "\t" + str(vog_name) + "\t" + str(annotations[n+11]) + "\t" + str(annotations[n+12]) + "\t" + str(annotations[n+13])
						anno_dict.update({annotations[n]:temp_result})
						n += 14
					for item in final:
						output.write(str(item).replace("$~&", " ").replace('^@%','"') + "\t" + str(anno_dict[item]) + "\n")
						gb_out.write(str(item).replace("$~&", " ").replace('^@%','"') + "\t" + str(genbank_dict_full[item]) + "\n")
						if item in AMG_dict:
							AMG_out.write(str(item).replace("$~&", " ").replace('^@%','"') + "\t" + str(AMG_dict[item]) + "\n")

with open(str(path)+'unmodified_VIBRANT_results.' + str(base) + '.txt', 'r') as results:
	with open(str(path)+'VIBRANT_results.' + str(base) + '.txt', 'w') as Vresults:
		results = results.read().replace('\n','\t').split('\t')
		if results[-1] == '':
			results = results[:-1]
		result_dict = {}
		prophage_genomes_frags = []
		prophage_genomes_base = []
		prophage_genomes = []
		phage_genomes = []

		n = 0
		while n < len(results):
			temp_result = str(results[n+1]).replace("$~&", " ").replace('^@%','"') + "\t" + str(results[n+2]) + "\t" + str(results[n+3]) + "\t" + str(results[n+4]) + "\t" + str(results[n+5]) + "\t" + str(results[n+6]) + "\t" + str(results[n+7]) + "\t" + str(results[n+8]) + "\t" + str(results[n+9]) + "\t" + str(results[n+10]) + "\t" + str(results[n+11]) + "\t" + str(results[n+12]) + "\t" + str(results[n+13]) + "\t" + str(results[n+14]) + "\t" + str(results[n+15]) + "\t" + str(results[n+16]) + "\t" + str(results[n+17]) + "\t" + str(results[n+18]) + "\t" + str(results[n+19]) + "\t" + str(results[n+20]) + "\t" + str(results[n+21]) + "\t" + str(results[n+22]) + "\t" + str(results[n+23]) + "\t" + str(results[n+24]) + "\t" + str(results[n+25]) + "\t" + str(results[n+26]) + "\t" + str(results[n+27])
			result_dict.update({results[n]:temp_result})
			n += 28
		for item in viral_genomes:
			Vresults.write(str(item).replace("$~&", " ").replace('^@%','"') + "\t" + str(result_dict[item]) + "\n")
			if "fragment" in str(item):
				prophage_genomes_frags.append(item)
				prophage_genomes_base.append(item.rsplit("_",2)[0])
			elif "fragment" not in str(item) and str(item) in prophage_integrase:
				prophage_genomes.append(item)
			elif "fragment" not in str(item) and str(item) not in prophage_integrase:
				phage_genomes.append(item)

prophage_genomes_frags = list(set(prophage_genomes_frags))
prophage_genomes_base = list(set(prophage_genomes_base))

######################## Write proteins, genes, genomes ########################
db_dict_protein = {}
db_dict_gene = {}
db_dict_protein_locs = {}
db_dict_protein_strand = {}
db_dict_gene_locs = {}
with open(in_proteins, 'r') as read_pass:
	with open(infile+'.phages_combined.faa', 'w') as combined_phages:
		with open(infile+'.phages_lysogenic.faa', 'w') as lysogenic_phages:
			with open(infile+'.phages_lytic.faa', 'w') as lytic_phages:
				for name, seq in SimpleFastaParser(read_pass):
					key = str(name.split(" # ")[0])
					first = name.split(" # ")[1]
					last = name.split(" # ")[2]
					strand = name.split(" # ")[3]
					db_dict_protein_locs.update({key:"("+str(first)+".."+str(last)+")"})
					db_dict_protein_strand.update({key:str(strand)})
					db_dict_protein.update({key:seq})
				for protein in final:
					protein_frag = str(protein.rsplit("_",3)[0]) + "_" + str(protein.rsplit("_",1)[1])
					if format == "prot":
						protein_replace = protein.replace("$~&", " ").replace('^@%','"')
						protein_frag_replace = protein_frag.replace("$~&", " ").replace('^@%','"')
						if "fragment" in protein:
							if protein_frag_replace in db_dict_protein:
								lysogenic_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein_frag_replace] + "\t" + db_dict_protein_strand[protein_frag_replace] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein_frag_replace] + "\n")
								combined_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein_frag_replace] + "\t" + db_dict_protein_strand[protein_frag_replace] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein_frag_replace] + "\n")
						if protein_replace in db_dict_protein and "fragment" not in protein:
							if str(protein.rsplit("_",1)[0]) not in prophage_integrase:
								lytic_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein_replace] + "\t" + db_dict_protein_strand[protein_replace] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein_replace] + "\n")
								combined_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein_replace] + "\t" + db_dict_protein_strand[protein_replace] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein_replace] + "\n")
							if str(protein.rsplit("_",1)[0]) in prophage_integrase:
								lysogenic_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein_replace] + "\t" + db_dict_protein_strand[protein_replace] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein_replace] + "\n")
								combined_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein_replace] + "\t" + db_dict_protein_strand[protein_replace] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein_replace] + "\n")
					if format == "nucl":
						if "fragment" in protein:
							if protein_frag in db_dict_protein:
								lysogenic_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein_frag] + "\t" + db_dict_protein_strand[protein_frag] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein_frag] + "\n")
								combined_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein_frag] + "\t" + db_dict_protein_strand[protein_frag] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein_frag] + "\n")
						if protein in db_dict_protein and "fragment" not in protein:
							if str(protein.rsplit("_",1)[0]) not in prophage_integrase:
								lytic_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein] + "\t" + db_dict_protein_strand[protein] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein] + "\n")
								combined_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein] + "\t" + db_dict_protein_strand[protein] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein] + "\n")
							if str(protein.rsplit("_",1)[0]) in prophage_integrase:
								lysogenic_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein] + "\t" + db_dict_protein_strand[protein] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein] + "\n")
								combined_phages.write(">" + protein.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_protein_locs[protein] + "\t" + db_dict_protein_strand[protein] + "\t" + genbank_dict[protein] + "\n" + db_dict_protein[protein] + "\n")

if format == "nucl":
	with open(infile+'.ffn', 'r') as read_pass:
		with open(infile+'.phages_combined.ffn', 'w') as combined_phages:
			with open(infile+'.phages_lysogenic.ffn', 'w') as lysogenic_phages:
				with open(infile+'.phages_lytic.ffn', 'w') as lytic_phages:
					for name, seq in SimpleFastaParser(read_pass):
						key = str(name.split(" # ")[0])
						first = name.split(" # ")[1]
						last = name.split(" # ")[2]
						db_dict_gene_locs.update({key:"("+str(first)+".."+str(last)+")"})
						db_dict_gene.update({key:seq})
					for gene in final:
						if "fragment" in gene:
							gene_frag = str(gene.rsplit("_",3)[0]) + "_" + str(gene.rsplit("_",1)[1])
							if gene_frag in db_dict_gene.keys():
								lysogenic_phages.write(">" + gene.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_gene_locs[gene_frag] + "\t" + genbank_dict[gene] + "\n" + db_dict_gene[gene_frag] + "\n")
								combined_phages.write(">" + gene.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_gene_locs[gene_frag] + "\t" + genbank_dict[gene] + "\n" + db_dict_gene[gene_frag] + "\n")
						if "fragment" not in gene:
							if gene in db_dict_gene.keys():
								if str(gene.rsplit("_",1)[0]) not in prophage_integrase:
									lytic_phages.write(">" + gene.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_gene_locs[gene] + "\t" + genbank_dict[gene] + "\n" + db_dict_gene[gene] + "\n")
									combined_phages.write(">" + gene.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_gene_locs[gene] + "\t" + genbank_dict[gene] + "\n" + db_dict_gene[gene] + "\n")
								if str(gene.rsplit("_",1)[0]) in prophage_integrase:
									lysogenic_phages.write(">" + gene.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_gene_locs[gene] + "\t" + genbank_dict[gene] + "\n" + db_dict_gene[gene] + "\n")
									combined_phages.write(">" + gene.replace("$~&", " ").replace('^@%','"') + "\t" + db_dict_gene_locs[gene] + "\t" + genbank_dict[gene] + "\n" + db_dict_gene[gene] + "\n")

	with open(input, 'r') as read_genomes:
		with open(infile+'.phages_combined.faa', 'r') as proteins:
			with open(infile+'.phages_combined.fna', 'w') as combined_phages:
				with open(infile+'.phages_lysogenic.fna', 'w') as lysogenic_phages:
					with open(infile+'.phages_lytic.fna', 'w') as lytic_phages:
						with open(in_proteins, 'r') as prodigal:
							definition_dict = {}
							remove_item = []
							for definition, temp in SimpleFastaParser(prodigal):
								key = str(definition.split(" # ")[0])
								first = definition.split(" # ")[1]
								last = definition.split(" # ")[2]
								definition_dict.update({key:str(first)+"~"+str(last)})
							for name, seq in SimpleFastaParser(read_genomes):
								if name in phage_genomes:
									lytic_phages.write(">" + str(name).replace("$~&", " ").replace('^@%','"') + "\n" + str(seq) + "\n")
									combined_phages.write(">" + str(name).replace("$~&", " ").replace('^@%','"') + "\n" + str(seq) + "\n")
								elif name in prophage_genomes:
									lysogenic_phages.write(">" + str(name).replace("$~&", " ").replace('^@%','"') + "\n" + str(seq) + "\n")
									combined_phages.write(">" + str(name).replace("$~&", " ").replace('^@%','"') + "\n" + str(seq) + "\n")
								elif name in prophage_genomes_base:
									for fragment in prophage_genomes_frags:
										if name == fragment.rsplit("_",2)[0]:
											frag_check = False
											genes = []
											for item in final:
												if item.rsplit("_",1)[0] == fragment:
													if item not in remove_item:
														genes.append(str(item.rsplit("_",3)[0])+"_"+str(item.rsplit("_",1)[1]))
														frag_check = True
														remove_item.append(item)
											if len(genes) > 1:
												gene_first = genes[0]
												gene_last = genes[-1]
												start = definition_dict[gene_first].split("~")[0]
												stop = definition_dict[gene_last].split("~")[1]
												seq_list = list(seq)
												sequence = seq_list[int(start)-1:int(stop)-1]
												sequence = "".join(sequence)
												lysogenic_phages.write(">" + str(fragment).replace("$~&", " ").replace('^@%','"') + "\n" + str(sequence) + "\n")
												combined_phages.write(">" + str(fragment).replace("$~&", " ").replace('^@%','"') + "\n" + str(sequence) + "\n")

################################### End of analysis ############################
####### remove files
subprocess.call('rm '+str(path)+'temp2_VIBRANT_annotations.' + str(base) + '.txt'+' 2>/dev/null', shell=True)
subprocess.call('rm '+str(path)+'temp_VIBRANT_results.' + str(base) + '.txt'+' 2>/dev/null', shell=True)
subprocess.call('rm '+str(path)+str(base)+'.pass.faa'+' 2>/dev/null', shell=True)
subprocess.call('rm '+str(path)+str(base)+'.appended*faa'+' 2>/dev/null', shell=True)
subprocess.call('rm '+str(path)+str(base)+'.master.txt'+' 2>/dev/null', shell=True)


                                                               ##')
                                                             ##  ##')
                                                           ##      ##')
######   ##  ##     ##     #######   ######    #####       ##      ##')
##  ##   ##  ##   ##  ##   ##        ##       ##             ##  ##')
######   ######   ######   ##  ###   ######    ###             ##')
##       ##  ##   ##  ##   ##   ##   ##           ##           ##')
##       ##  ##   ##  ##   #######   ######   #####            ##')
                                                            #  ##  #')
                                                           # # ## # #')
                                                          #   #  #   #')
                                                         #            #')

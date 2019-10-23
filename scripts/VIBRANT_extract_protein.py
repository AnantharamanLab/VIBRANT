#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison, 2019

# VIBRANT v1.0.1
# Virus Identification By iteRative ANnoTation

# Usage: see VIBRANT_run.py

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

with open(str(sys.argv[1]), 'r') as accnos:
    with open(str(sys.argv[2]), 'w') as fasta:
        with open(str(sys.argv[3]), 'r') as db:
            accnos_list = accnos.read().split("\n")
            if accnos_list[-1] == '':
                accnos_list = accnos_list[:-1]
            for name, seq in SimpleFastaParser(db):
                temp = name.split(" # ",1)[0]
                if temp.rsplit("_",1)[0] in accnos_list:
                    fasta.write(">" + name + "\n" + seq + "\n")

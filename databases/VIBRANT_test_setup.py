#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison, 2019

# VIBRANT
# Virus Identification By iteRative ANnoTation

# Usage: $ python3 VIBRANT_test_setup.py

try:
    import sys
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from datetime import date
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
    from collections import OrderedDict
    import numpy as np
    from sklearn import preprocessing
    from sklearn.neural_network import MLPClassifier
    from sklearn.metrics import classification_report,confusion_matrix
    import pickle
    import os
    import sklearn
except Exception as e:
    print("\n")
    print("VIBRANT Error: Please install VIBRANT dependancy module(s). See below for the encountered error:")
    print("\n")
    print(e)
    print("\n")
    exit()

error = 0
if str(sklearn.__version__) != '0.21.3':
    if int(str(sklearn.__version__).split(".")[1]) > 21 or int(str(sklearn.__version__).split(".")[2]) > 3:
        print()
        print('VIBRANT Caution: running a version of Scikit-Learn higher than v0.21.3 may cause issues. With pip you can update by running "pip install --upgrade scikit-learn==0.21.3".')
        print()
        error += 1
    elif int(str(sklearn.__version__).split(".")[1]) < 21 or int(str(sklearn.__version__).split(".")[2]) < 3:
        print()
        print('VIBRANT Caution: running a version of Scikit-Learn lower than v0.21.3 will likely cause issues. With pip you can update by running "pip install --upgrade scikit-learn==0.21.3".')
        print()
        error += 1


if str(np.version.version) != '1.17.0':
    if int(str(np.version.version).split(".")[0]) < 1 or int(str(np.version.version).split(".")[1]) < 17:
        print()
        print('VIBRANT Caution: running a version of Numpy lower than v1.17.0 will likely cause issues. With pip you can update by running "pip install --upgrade numpy==1.17.0".')
        print()
        error += 1

parent_path = str(os.path.dirname(os.path.abspath(__file__)).rsplit("/",1)[0])

listing = 'grep -c "NAME" ' + str(parent_path) + '/databases/KEGG_profiles_prokaryotes.HMM'
listing_shell = subprocess.check_output(listing, shell=True)
listing = str(listing_shell.strip()).split("'")[1]
kegg = str(listing.split("\\n")[0].split(" ",1)[0])

listing = 'grep -c "NAME" ' + str(parent_path) + '/databases/Pfam-A_phage_v32.HMM'
listing_shell = subprocess.check_output(listing, shell=True)
listing = str(listing_shell.strip()).split("'")[1]
pfam_phage = str(listing.split("\\n")[0].split(" ",1)[0])

listing = 'grep -c "NAME" ' + str(parent_path) + '/databases/Pfam-A_plasmid_v32.HMM'
listing_shell = subprocess.check_output(listing, shell=True)
listing = str(listing_shell.strip()).split("'")[1]
pfam_plasmid = str(listing.split("\\n")[0].split(" ",1)[0])

listing = 'grep -c "NAME" ' + str(parent_path) + '/databases/Pfam-A_v32.HMM	'
listing_shell = subprocess.check_output(listing, shell=True)
listing = str(listing_shell.strip()).split("'")[1]
pfam = str(listing.split("\\n")[0].split(" ",1)[0])

listing = 'grep -c "NAME" ' + str(parent_path) + '/databases/VOGDB94_phage.HMM'
listing_shell = subprocess.check_output(listing, shell=True)
listing = str(listing_shell.strip()).split("'")[1]
vog = str(listing.split("\\n")[0].split(" ",1)[0])

try:
    subprocess.check_output("prodigal -h 2>/dev/null", shell=True)
except Exception as e:
    print("VIBRANT Error: Prodigal is not installed or not in $PATH. Please install Prodigal or add to $PATH.")
    print("\n")
    exit()

try:
    subprocess.check_output("hmmsearch -h 2>/dev/null", shell=True)
except Exception as e:
    print("VIBRANT Error: HMMER3 is not installed or not in $PATH. Please install HMMER3 or add to $PATH.")
    print("\n")
    exit()


if str(kegg) != "10033":
    print("VIBRANT Error: it looks like the KEGG HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1
if str(pfam_phage) != "894":
    print("VIBRANT Error: it looks like the Pfam phage HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1
if str(pfam_plasmid) != "202":
    print("VIBRANT Error: it looks like the Pfam plasmid HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1
if str(pfam) != "17929":
    print("VIBRANT Error: it looks like the Pfam HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1
if str(vog) != "19182":
    print("VIBRANT Error: it looks like the VOG HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1

if error > 0:
    exit()

if error == 0:
    print("\n")
    print("VIBRANT v1.0.0 is good to go!")
    print("\n")

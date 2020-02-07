#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison

# VIBRANT v1.1.0
# Virus Identification By iteRative ANnoTation
# Release date: Feb 7 2020

# Usage: $ python3 VIBRANT_setup.py

import subprocess
from subprocess import PIPE
import time
import os
import argparse
import sys
from io import StringIO
import logging

vibrant = argparse.ArgumentParser(description='Usage: python3 VIBRANT_setup.py.')
vibrant.add_argument('--version', action='version', version='VIBRANT v1.1.0')
vibrant.add_argument('-test', action='store_true', help='use this setting if you only want to test downloads/dependencies [default=off]')
vibrant.add_argument('-force', action='store_true', help='use this setting if this script is exiting and telling you to install a package/program that you know is already installed [default=off]')
args = vibrant.parse_args()

parent_path = str(os.path.dirname(os.path.abspath(__file__)).rsplit("/",1)[0])

if args.test == False:
    print()
    print("This script will download, extract subsets and press HMM profiles for VIBRANT.")
    print("This process will require 20GB of temporary free storage space, but the final size requirement is ~11.2GB in the form of pressed HMM databases.")
    print("Please be patient. This only needs to be run once and will take a few minutes.")
    print("Logger started. Check log file for messages and errors.")

    logging.basicConfig(filename=str(parent_path) + '/databases/VIBRANT_setup.log', level=logging.INFO, format='%(message)s')
    logging.info = logging.info
    logging.info("This script will download, extract subsets and press HMM profiles for VIBRANT.")
    logging.info("This process will require 20GB of temporary free storage space, but the final size requirement is ~11.2GB in the form of pressed HMM databases.")
    logging.info("Please be patient. This only needs to be run once and will take a few minutes.")
    logging.info('')

    if args.force == False:
        try:
            subprocess.check_output("wget --help 2>/dev/null", shell=True)
        except Exception as e:
            logging.info("VIBRANT Error: wget is not installed. Please install wget.")
            logging.info("\n")
            exit()

        try:
            subprocess.check_output("hmmfetch -h 2>/dev/null", shell=True)
        except Exception as e:
            logging.info("VIBRANT Error: HMMER3 is not installed or not in $PATH. Please install HMMER3.")
            logging.info("\n")
            exit()

        try:
            subprocess.check_output("gunzip -h 2>/dev/null", shell=True)
        except Exception as e:
            logging.info("VIBRANT Error: gunzip is not installed. Please install gunzip.")
            logging.info("\n")
            exit()

        try:
            subprocess.check_output("tar --help 2>/dev/null", shell=True)
        except Exception as e:
            logging.info("VIBRANT Error: tar is not installed. Please install tar.")
            logging.info("\n")
            exit()

    logging.info("Verifying Pfam, KEGG and VOG source websites are available for download ...")
    c1 = subprocess.Popen("wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz --spider -nv --no-hsts -t 1", shell=True, stderr=PIPE)
    c2 = subprocess.Popen("wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz --spider -nv --no-hsts -t 1", shell=True, stderr=PIPE)
    c3 = subprocess.Popen("wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz --spider -nv --no-hsts -t 1", shell=True, stderr=PIPE)
    vog_out = c1.stderr.read()
    pfam_out = c2.stderr.read()
    kegg_out = c3.stderr.read()
    if "failed" in str(vog_out):
        logging.info("\n")
        logging.info("Error: looks like the VOG source website is not working. If the link is not working you may need to try later. Sorry, this is out of the control of VIBRANT.")
        logging.info("Here is the source link: http://fileshare.csb.univie.ac.at/vog/vog94/")
        logging.info("\n")
        exit()
    if "failed" in str(pfam_out):
        logging.info("\n")
        logging.info("Error: looks like the Pfam source website is not working. If the link is not working you may need to try later. Sorry, this is out of the control of VIBRANT.")
        logging.info("Here is the source link: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/")
        logging.info("\n")
        exit()
    if "failed" in str(kegg_out):
        logging.info("\n")
        logging.info("Error: looks like the KEGG source website is not working. If the link is not working you may need to try later. Sorry, this is out of the control of VIBRANT.")
        logging.info("Here is the source link: ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/")
        logging.info("\n")
        exit()


    logging.info('')
    logging.info("Downloading HMM profiles for Pfam, KEGG and VOG from their source websites ...")
    d1 = subprocess.Popen('wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz -q --no-hsts', shell=True)
    d2 = subprocess.Popen('wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz -q --no-hsts', shell=True)
    d3 = subprocess.Popen('wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz -q --no-hsts', shell=True)
    d1.wait()
    d2.wait()
    d3.wait()
    time.sleep(0.1)

    logging.info('')
    logging.info("Unzipping profiles ...")
    u1 = subprocess.Popen('tar -xzf vog.hmm.tar.gz', shell=True)
    u2 = subprocess.Popen('gunzip Pfam-A.hmm.gz', shell=True)
    u3 = subprocess.Popen('tar -xzf profiles.tar.gz', shell=True)
    u1.wait()
    u2.wait()
    u3.wait()
    time.sleep(0.1)

    logging.info('')
    logging.info("Concatenating individual profiles ...")
    c1 = subprocess.Popen('for v in VOG*.hmm; do cat $v >> vog_temp.HMM; done', shell=True)
    c2 = subprocess.Popen('for k in profiles/K*.hmm; do cat $k >> kegg_temp.HMM; done', shell=True)
    c1.wait()
    c2.wait()
    time.sleep(0.1)
    c3 = subprocess.Popen('rm VOG0*.hmm', shell=True)
    c4 = subprocess.Popen('rm VOG1*.hmm', shell=True)
    c5 = subprocess.Popen('rm VOG2*.hmm', shell=True)
    c6 = subprocess.Popen('rm -R profiles', shell=True)
    c3.wait()
    c4.wait()
    c5.wait()
    c6.wait()
    time.sleep(0.1)

    logging.info('')
    logging.info("Extracting profiles used for VIBRANT ...")
    e1 = subprocess.Popen('hmmfetch -o VOGDB94_phage.HMM -f vog_temp.HMM profile_names/VIBRANT_vog_profiles.txt >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    e2 = subprocess.Popen('hmmfetch -o KEGG_profiles_prokaryotes.HMM -f kegg_temp.HMM profile_names/VIBRANT_kegg_profiles.txt >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    e3 = subprocess.Popen('hmmfetch -o Pfam-A_plasmid_v32.HMM -f Pfam-A.hmm profile_names/VIBRANT_pfam-plasmid_profiles.txt >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    e4 = subprocess.Popen('hmmfetch -o Pfam-A_phage_v32.HMM -f Pfam-A.hmm profile_names/VIBRANT_pfam-phage_profiles.txt >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    e1.wait()
    e2.wait()
    e3.wait()
    e4.wait()
    time.sleep(0.1)
    r1 = subprocess.Popen('rm vog_temp.HMM kegg_temp.HMM vog.hmm.tar.gz profiles.tar.gz', shell=True)
    r2 = subprocess.Popen('mv Pfam-A.hmm Pfam-A_v32.HMM', shell=True)
    r1.wait()
    r2.wait()
    time.sleep(0.1)

    logging.info('')
    logging.info("Pressing profiles used for VIBRANT ...")
    p1 = subprocess.Popen('hmmpress VOGDB94_phage.HMM >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    p2 = subprocess.Popen('hmmpress KEGG_profiles_prokaryotes.HMM >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    p3 = subprocess.Popen('hmmpress Pfam-A_plasmid_v32.HMM >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    p4 = subprocess.Popen('hmmpress Pfam-A_phage_v32.HMM >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    p5 = subprocess.Popen('hmmpress Pfam-A_v32.HMM >> '+ str(parent_path) + '/databases/VIBRANT_setup.log', shell=True)
    p1.wait()
    p2.wait()
    p3.wait()
    p4.wait()
    p5.wait()
    time.sleep(0.1)

    logging.info('')
    logging.info("Done with databases. Several new databases are now in this folder.")
    logging.info('')
    logging.info("Verying correct dependency versions ...")
    logging.info('')
    #
    #
    #
if args.test == True:
    print()
    print("Verying correct dependency versions ...")
    print("Logger started. Check log file for messages and errors.")
    logging.basicConfig(filename=str(parent_path) + '/databases/VIBRANT_test_setup.log', level=logging.INFO, format='%(message)s')
    logging.info = logging.info

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
    logging.info("\n")
    logging.info("VIBRANT Error: Please install VIBRANT dependancy module(s). See below for the encountered error:")
    logging.info('')
    logging.info(e)
    logging.info('')
    exit()

error = 0
if str(sklearn.__version__) != '0.21.3':
    if int(str(sklearn.__version__).split(".")[1]) > 21 or int(str(sklearn.__version__).split(".")[2]) > 3:
        logging.info('')
        logging.info('VIBRANT Caution: running a version of Scikit-Learn higher than v0.21.3 may cause issues. With pip you can update by running "pip install --upgrade scikit-learn==0.21.3".')
        logging.info('')
        error += 1
    elif int(str(sklearn.__version__).split(".")[1]) < 21 or int(str(sklearn.__version__).split(".")[2]) < 3:
        logging.info('')
        logging.info('VIBRANT Caution: running a version of Scikit-Learn lower than v0.21.3 will likely cause issues. With pip you can update by running "pip install --upgrade scikit-learn==0.21.3".')
        logging.info('')
        error += 1
if str(np.version.version) != '1.17.0':
    if int(str(np.version.version).split(".")[0]) < 1 or int(str(np.version.version).split(".")[1]) < 17:
        logging.info('')
        logging.info('VIBRANT Caution: running a version of Numpy lower than v1.17.0 will likely cause issues. With pip you can update by running "pip install --upgrade numpy==1.17.0".')
        logging.info('')
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
    logging.info("VIBRANT Error: Prodigal is not installed or not in $PATH. Please install Prodigal or add to $PATH.")
    logging.info("\n")
    exit()

try:
    subprocess.check_output("hmmsearch -h 2>/dev/null", shell=True)
except Exception as e:
    logging.info("VIBRANT Error: HMMER3 is not installed or not in $PATH. Please install HMMER3 or add to $PATH.")
    logging.info("\n")
    exit()

if str(kegg) != "10033":
    logging.info("VIBRANT Error: it looks like the KEGG HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1
if str(pfam_phage) != "894":
    logging.info("VIBRANT Error: it looks like the Pfam phage HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1
if str(pfam_plasmid) != "202":
    logging.info("VIBRANT Error: it looks like the Pfam plasmid HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1
if str(pfam) != "17929":
    logging.info("VIBRANT Error: it looks like the Pfam HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1
if str(vog) != "19182":
    logging.info("VIBRANT Error: it looks like the VOG HMM profiles were not downloaded/pressed correctly. Try re-running the VIBRANT_setup.py script.")
    error += 1

if error > 0:
    exit()

if error == 0:
    logging.info('')
    logging.info("VIBRANT v1.1.0 is good to go!")
    logging.info('See example_data/ for quick test files.')
    logging.info('')

print()
print("VIBRANT v1.1.0 is good to go!")
print("See example_data/ for quick test files.")
print()

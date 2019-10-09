#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison, 2019

# VIBRANT
# Virus Identification By iteRative ANnoTation

# Usage: $ python3 VIBRANT_setup.py

import subprocess
from subprocess import PIPE
import time
import os
import argparse
import sys
from io import StringIO

vibrant = argparse.ArgumentParser(description='Usage: python3 VIBRANT_setup.py.')
vibrant.add_argument('-force', action='store_true', help='use this setting if this script is exiting and telling you to install a package/program that you know is already installed [default=off]')
args = vibrant.parse_args()

parent_path = str(os.path.dirname(os.path.abspath(__file__)).rsplit("/",1)[0])

print("\n")
print("This script will download, extract subsets and press HMM profiles for VIBRANT.")
print("This process will require 20GB of temporary free storage space, but the final size requirement is ~11.2GB in the form of pressed HMM databases.")
print("Please be patient. This only needs to be run once and will take a few minutes.")
print("\n")

if args.force == False:
    try:
        subprocess.check_output("wget --help 2>/dev/null", shell=True)
    except Exception as e:
        print("VIBRANT Error: wget is not installed. Please install wget.")
        print("\n")
        exit()

    try:
        subprocess.check_output("hmmfetch -h 2>/dev/null", shell=True)
    except Exception as e:
        print("VIBRANT Error: HMMER3 is not installed or not in $PATH. Please install HMMER3.")
        print("\n")
        exit()

    try:
        subprocess.check_output("gunzip -h 2>/dev/null", shell=True)
    except Exception as e:
        print("VIBRANT Error: gunzip is not installed. Please install gunzip.")
        print("\n")
        exit()

    try:
        subprocess.check_output("tar --help 2>/dev/null", shell=True)
    except Exception as e:
        print("VIBRANT Error: tar is not installed. Please install tar.")
        print("\n")
        exit()

print("Verifying Pfam, KEGG and VOG source websites are available for download ...")
c1 = subprocess.Popen("wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz --spider -nv --no-hsts -t 1", shell=True, stderr=PIPE)
c2 = subprocess.Popen("wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz --spider -nv --no-hsts -t 1", shell=True, stderr=PIPE)
c3 = subprocess.Popen("wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz --spider -nv --no-hsts -t 1", shell=True, stderr=PIPE)
vog_out = c1.stderr.read()
pfam_out = c2.stderr.read()
kegg_out = c3.stderr.read()
if "failed" in str(vog_out):
    print("\n")
    print("Error: looks like the VOG source website is not working. If the link is not working you may need to try later. Sorry, this is out of the control of VIBRANT.")
    print("Here is the source link: http://fileshare.csb.univie.ac.at/vog/vog94/")
    print("\n")
    exit()
if "failed" in str(pfam_out):
    print("\n")
    print("Error: looks like the Pfam source website is not working. If the link is not working you may need to try later. Sorry, this is out of the control of VIBRANT.")
    print("Here is the source link: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/")
    print("\n")
    exit()
if "failed" in str(kegg_out):
    print("\n")
    print("Error: looks like the KEGG source website is not working. If the link is not working you may need to try later. Sorry, this is out of the control of VIBRANT.")
    print("Here is the source link: ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/")
    print("\n")
    exit()


print("\n")
print("Downloading HMM profiles for Pfam, KEGG and VOG from their source websites ...")
d1 = subprocess.Popen('wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz -nv --no-hsts', shell=True)
d2 = subprocess.Popen('wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz -nv --no-hsts', shell=True)
d3 = subprocess.Popen('wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz -nv --no-hsts', shell=True)
d1.wait()
d2.wait()
d3.wait()
time.sleep(0.1)

print("\n")
print("Unzipping profiles ...")
u1 = subprocess.Popen('tar -xzf vog.hmm.tar.gz', shell=True)
u2 = subprocess.Popen('gunzip Pfam-A.hmm.gz', shell=True)
u3 = subprocess.Popen('tar -xzf profiles.tar.gz', shell=True)
u1.wait()
u2.wait()
u3.wait()
time.sleep(0.1)

print("\n")
print("Concatenating individual profiles ...")
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

print("\n")
print("Extracting profiles used for VIBRANT ...")
e1 = subprocess.Popen('hmmfetch -o VOGDB94_phage.HMM -f vog_temp.HMM profile_names/VIBRANT_vog_profiles.txt', shell=True)
e2 = subprocess.Popen('hmmfetch -o KEGG_profiles_prokaryotes.HMM -f kegg_temp.HMM profile_names/VIBRANT_kegg_profiles.txt', shell=True)
e3 = subprocess.Popen('hmmfetch -o Pfam-A_plasmid_v32.HMM -f Pfam-A.hmm profile_names/VIBRANT_pfam-plasmid_profiles.txt', shell=True)
e4 = subprocess.Popen('hmmfetch -o Pfam-A_phage_v32.HMM -f Pfam-A.hmm profile_names/VIBRANT_pfam-phage_profiles.txt', shell=True)
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

print("\n")
print("Pressing profiles used for VIBRANT ...")
p1 = subprocess.Popen('hmmpress VOGDB94_phage.HMM', shell=True)
p2 = subprocess.Popen('hmmpress KEGG_profiles_prokaryotes.HMM', shell=True)
p3 = subprocess.Popen('hmmpress Pfam-A_plasmid_v32.HMM', shell=True)
p4 = subprocess.Popen('hmmpress Pfam-A_phage_v32.HMM', shell=True)
p5 = subprocess.Popen('hmmpress Pfam-A_v32.HMM', shell=True)
p1.wait()
p2.wait()
p3.wait()
p4.wait()
p5.wait()
time.sleep(0.1)

print("\n")
print("Done. Several new databases are now in this folder.")
print("\n")
print("VIBRANT should be ready to go. You can verify this by running VIBRANT_test_setup.py within this folder (databases/)")
print("\n")
#
#
#

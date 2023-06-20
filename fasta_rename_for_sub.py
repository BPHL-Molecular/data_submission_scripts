#!/usr/bin/env python

#Author: Sarah Schmedes
#Email: sarah.schmedes@flhealth.gov

'''
This script will take in fasta assemblies, rename the headers and output a concatenated multi-fasta.
This is compatible with submissions to NCBI GenBank and GISAID.
'''

import os
import subprocess
import sys
import datetime
import argparse
import pandas as pd
from Bio import SeqIO

#Parse arguments, get path of renaming file
parser = argparse.ArgumentParser(usage='fasta_rename_for_sub.py --assem_dir <assemblies dir> --names <name file> ')
parser.add_argument('--assem_dir', help='path to assembly directory as input', required=True)
parser.add_argument('--names', help='path to sample renaming file', required=True)

if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

assem_dir = args.assem_dir
names = args.names
cwd = os.getcwd() + '/'

if assem_dir[-1] != '/':
    assem_dir = assem_dir + '/'

output_dir = cwd + 'fastas_for_submission/'
subprocess.run('mkdir -p ' + output_dir, shell=True, check=True)

#Read in names file
names_df = pd.read_table(names, sep ="\t", header=None)
lab_name = list(names_df[0])
public_name = list(names_df[1])
name_key = dict(zip(lab_name, public_name))

#Make temp fasta dir
subprocess.run('mkdir temp_fastas', shell=True, check=True)

#Find assembly in assembly directory and prepare for submission
for n in lab_name:
    #Get fasta file name
    proc_f = subprocess.run('ls ' + assem_dir + n + '*.fa*', shell=True, capture_output=True, text=True, check=True)
    fasta_file = proc_f.stdout.rstrip()
    #Read in each fasta using biopython and print back to file to get one line of sequence.
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        sn = seq_id.split("|")
        sn = sn[0]
        assert n == sn, 'File sample name does not match fasta header file name. Sample:' + n + '; header:' + sn
        seq = str(record.seq)
        with open('temp_fastas/' + n + '.fasta', 'w') as temp_fasta:
            temp_fasta.write('>' + seq_id + "\n")
            temp_fasta.write(seq)
    #Write header to new fasta and then write sequence (removing leading Ns and fold every 75 characters)
    subprocess.run('echo ">' + name_key[n] + '" > ' + output_dir + n + '.fasta', shell=True, check=True)
    subprocess.run('grep -v ">" temp_fastas/' + n + '.fasta | sed \'s/^N*N//g\' | fold -w 75 >> ' + output_dir + n + '.fasta', shell=True, check=True)

#Concatenate fastas for each submission file
subprocess.run('cat ' + output_dir + '*.fasta > ' + output_dir + datetime.date.today().strftime('%Y%m%d') + '_all.fasta', shell=True, check=True)

#Remove separated assemblies
subprocess.run('rm -r temp_fastas/', shell=True, check=True)

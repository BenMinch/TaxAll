import subprocess
import csv
import pandas as pd
import os, re, sys, argparse

parser = argparse.ArgumentParser(description='''Assign Taxonomy to Sequences''')
parser.add_argument('-i', '--input', help='''Input fasta PolB Seqs Folder''', required=True)
args = parser.parse_args()

input = args.input

for file in os.listdir(input):
    input_file = os.path.join(input, file)
    output_file = re.sub('.fa', '_tax', input_file)
    cmd='python TaxAll.py -i '+input_file+' -o '+output_file+ ' -d custom_refseq_ncldv2.dmnd -names total_names.dmp -nodes all_nodes.dmp'
    subprocess.call(cmd, shell=True)

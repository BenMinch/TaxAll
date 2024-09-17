import subprocess
import csv
import pandas as pd
import os, re, sys, argparse

parser = argparse.ArgumentParser(description='''Assign Taxonomy to Sequences''')
parser.add_argument('-i', '--input', help='''Input fasta PolB Seqs''', required=True)
parser.add_argument('-o', '--output', help='''Output Name ''', required=True)
parser.add_argument('-d', '--database', help='''Database Name''', required=True)
parser.add_argument('-names', '--names', help='''Names File''', required=True)
parser.add_argument('-nodes', '--nodes', help='''Nodes File''', required=True)
parser.add_argument('-n', '--n', help='''Number of Hits you want to see''', required=False, default='3')

args = parser.parse_args()

input = args.input
output = args.output
database = args.database
names = args.names
nodes = args.nodes
n = args.n



def run_diamond_blastp(query_file, database_file, output_file,n):
    cmd = f"diamond blastp -q {query_file} -d {database_file} -o {output_file} --outfmt 6 qseqid sseqid staxids pident --max-target-seqs {n}"
    subprocess.run(cmd, shell=True, check=True)

# Replace with your input FASTA file, custom database file, and output CSV file
fasta_file = input
database_file = database
output_file = output+'blast_results.txt'
output_csv = output+'.csv'

# Run DIAMOND BLAST
run_diamond_blastp(fasta_file, database_file, output_file,n)


data=pd.read_csv(nodes,sep='\t',header=None)

#only keep column 0, 2, and 4
names= pd.read_csv(names,sep='\t',header=None)

data= data[[0,2,4]]
data.columns= ['taxid','parent_taxid','rank']
#get all unique ranks


names= names[[0,2,6]]
names.columns= ['taxid','name','type']

names_sci= names[names['type']=='scientific name']
#Get species name
blast_data=pd.read_csv(output_file,sep='\t',header=None)

blast_data.columns= ['query','subject','taxid','pident']

#Get Scientific Name
blast_data['species']=''
blast_data['genus']=''
blast_data['family']=''
blast_data['order']=''
blast_data['class']=''
blast_data['phylum']=''
blast_data['kingdom']=''
blast_data['superkingdom']=''
blast_data['clade']=''
for i in range(len(blast_data)):
    kingdom_value=''
    class_value=''
    order_value=''
    family_value=''
    genus_value=''
    phylum_value=''
    superkingdom_value=''
    next_value=''
    next2_value=''
    next3_value=''
    next4_value=''
    next5_value=''
    next6_value=''
    next7_value=''
    next8_value=''
    next9_value=''
    next10_value=''
    next11_value=''
    next12_value=''
    next13_value=''
    next14_value=''
    next15_value=''
    next16_value=''
    next17_value=''
    next18_value=''
    next19_value=''
    next20_value=''
    taxid= blast_data['taxid'][i]
    print(taxid)
    #if taxid==2: skip
    if taxid==2:
        continue
    if taxid in names_sci['taxid'].values:
        species= names_sci[taxid==names_sci['taxid']]['name'].values[0]
        species_value= data[taxid==data['taxid']]['rank'].values[0]
        if species_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
            blast_data[species_value][i]= species
        else:
            pass
        genus= data[taxid==data['taxid']]['parent_taxid'].values[0]
        genus_name= names_sci[genus==names_sci['taxid']]['name'].values[0]
        genus_value= data[genus==data['taxid']]['rank'].values[0]
        if genus_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
            blast_data[genus_value][i]= genus_name
        else:
            pass
        if genus_value in ['superkingdom','']:
            pass
        else:
            family= data[genus==data['taxid']]['parent_taxid'].values[0]
            family_name= names_sci[family==names_sci['taxid']]['name'].values[0]
            family_value= data[family==data['taxid']]['rank'].values[0]
            if family_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[family_value][i]= family_name
            else:
                pass
        if family_value in ['superkingdom','']:
            pass
        else:
            order= data[family==data['taxid']]['parent_taxid'].values[0]
            order_name= names_sci[order==names_sci['taxid']]['name'].values[0]
            order_value= data[order==data['taxid']]['rank'].values[0]
            if order_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[order_value][i]= order_name
            else:
                pass
        if order_value in ['superkingdom','']:
            pass
        else:
            class_= data[order==data['taxid']]['parent_taxid'].values[0]
            class_name= names_sci[class_==names_sci['taxid']]['name'].values[0]
            class_value= data[class_==data['taxid']]['rank'].values[0]
            if class_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[class_value][i]= class_name
            else:
                pass
        if class_value in ['superkingdom','']:
            pass
        else:
            phylum= data[class_==data['taxid']]['parent_taxid'].values[0]
            phylum_name= names_sci[phylum==names_sci['taxid']]['name'].values[0]
            phylum_value= data[phylum==data['taxid']]['rank'].values[0]
            if phylum_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[phylum_value][i]= phylum_name
            else:
                pass
        if phylum_value in ['superkingdom','']:
            pass
        else:
            kingdom= data[phylum==data['taxid']]['parent_taxid'].values[0]
            kingdom_name= names_sci[kingdom==names_sci['taxid']]['name'].values[0]
            kingdom_value= data[kingdom==data['taxid']]['rank'].values[0]
            if kingdom_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[kingdom_value][i]= kingdom_name
            else:
                pass
        if kingdom_value in ['superkingdom','']:
            pass
        else:
            superkingdom= data[kingdom==data['taxid']]['parent_taxid'].values[0]
            superkingdom_name= names_sci[superkingdom==names_sci['taxid']]['name'].values[0]
            superkingdom_value= data[superkingdom==data['taxid']]['rank'].values[0]
            if superkingdom_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[superkingdom_value][i]= superkingdom_name
            else:
                pass
        if superkingdom_value in ['superkingdom','']:
            pass
        else:
            next_= data[superkingdom==data['taxid']]['parent_taxid'].values[0]
            next_name= names_sci[next_==names_sci['taxid']]['name'].values[0]
            next_value= data[next_==data['taxid']]['rank'].values[0]
            if next_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next_value][i]= next_name
            else:
                pass
        if next_value in ['superkingdom','']:
            pass
        else:
            next2_= data[next_==data['taxid']]['parent_taxid'].values[0]
            next2_name= names_sci[next2_==names_sci['taxid']]['name'].values[0]
            next2_value= data[next2_==data['taxid']]['rank'].values[0]
            if next2_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next2_value][i]= next2_name
            else:
                pass
        if next2_value in ['superkingdom','']:
            pass
        else:
            next3_= data[next2_==data['taxid']]['parent_taxid'].values[0]
            next3_name= names_sci[next3_==names_sci['taxid']]['name'].values[0]
            next3_value= data[next3_==data['taxid']]['rank'].values[0]
            if next3_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next3_value][i]= next3_name
            else:
                pass
        if next3_value in ['superkingdom','']:
            pass
        else:
            next4_= data[next3_==data['taxid']]['parent_taxid'].values[0]
            next4_name= names_sci[next4_==names_sci['taxid']]['name'].values[0]
            next4_value= data[next4_==data['taxid']]['rank'].values[0]
            if next4_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next4_value][i]= next4_name
            else:
                pass
        if next4_value in ['superkingdom','']:
            pass
        else:
            next5_= data[next4_==data['taxid']]['parent_taxid'].values[0]
            next5_name= names_sci[next5_==names_sci['taxid']]['name'].values[0]
            next5_value= data[next5_==data['taxid']]['rank'].values[0]
            if next5_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next5_value][i]= next5_name
            else:
                pass
        if next5_value in ['superkingdom','']:
            pass
        else:
            next6_= data[next5_==data['taxid']]['parent_taxid'].values[0]
            next6_name= names_sci[next6_==names_sci['taxid']]['name'].values[0]
            next6_value= data[next6_==data['taxid']]['rank'].values[0]
            if next6_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next6_value][i]= next6_name
            else:
                pass
        if next6_value in ['superkingdom','']:
            pass
        else:
            next7_= data[next6_==data['taxid']]['parent_taxid'].values[0]
            next7_name= names_sci[next7_==names_sci['taxid']]['name'].values[0]
            next7_value= data[next7_==data['taxid']]['rank'].values[0]
            if next7_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next7_value][i]= next7_name
            else:
                pass
        if next7_value in ['superkingdom','']:
            pass
        else:
            next8_= data[next7_==data['taxid']]['parent_taxid'].values[0]
            next8_name= names_sci[next8_==names_sci['taxid']]['name'].values[0]
            next8_value= data[next8_==data['taxid']]['rank'].values[0]
            if next8_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next8_value][i]= next8_name
            else:
                pass
        if next8_value in ['superkingdom','']:
            pass
        else:
            next9_= data[next8_==data['taxid']]['parent_taxid'].values[0]
            next9_name= names_sci[next9_==names_sci['taxid']]['name'].values[0]
            next9_value= data[next9_==data['taxid']]['rank'].values[0]
            if next9_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next9_value][i]= next9_name
            else:
                pass
        if next9_value in ['superkingdom','']:
            pass
        else:
            next10_= data[next9_==data['taxid']]['parent_taxid'].values[0]
            next10_name= names_sci[next10_==names_sci['taxid']]['name'].values[0]
            next10_value= data[next10_==data['taxid']]['rank'].values[0]
            if next10_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next10_value][i]= next10_name
            else:
                pass
        if next10_value in ['superkingdom','']:
            pass
        else:
            next11_= data[next10_==data['taxid']]['parent_taxid'].values[0]
            next11_name= names_sci[next11_==names_sci['taxid']]['name'].values[0]
            next11_value= data[next11_==data['taxid']]['rank'].values[0]
            if next11_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next11_value][i]= next11_name
            else:
                pass
        if next11_value in ['superkingdom','']:
            pass
        else:
            next12_= data[next11_==data['taxid']]['parent_taxid'].values[0]
            next12_name= names_sci[next12_==names_sci['taxid']]['name'].values[0]
            next12_value= data[next12_==data['taxid']]['rank'].values[0]
            if next12_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next12_value][i]= next12_name
            else:
                pass
        if next12_value in ['superkingdom','']:
            pass
        else:
            next13_= data[next12_==data['taxid']]['parent_taxid'].values[0]
            next13_name= names_sci[next13_==names_sci['taxid']]['name'].values[0]
            next13_value= data[next13_==data['taxid']]['rank'].values[0]
            if next13_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next13_value][i]= next13_name
            else:
                pass
        if next13_value in ['superkingdom','']:
            pass
        else:
            next14_= data[next13_==data['taxid']]['parent_taxid'].values[0]
            next14_name= names_sci[next14_==names_sci['taxid']]['name'].values[0]
            next14_value= data[next14_==data['taxid']]['rank'].values[0]
            if next14_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next14_value][i]= next14_name
            else:
                pass
        if next14_value in ['superkingdom','']:
            pass
        else:
            next15_= data[next14_==data['taxid']]['parent_taxid'].values[0]
            next15_name= names_sci[next15_==names_sci['taxid']]['name'].values[0]
            next15_value= data[next15_==data['taxid']]['rank'].values[0]
            if next15_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next15_value][i]= next15_name
            else:
                pass
        if next15_value in ['superkingdom','']:
            pass
        else:
            next16_= data[next15_==data['taxid']]['parent_taxid'].values[0]
            next16_name= names_sci[next16_==names_sci['taxid']]['name'].values[0]
            next16_value= data[next16_==data['taxid']]['rank'].values[0]
            if next16_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next16_value][i]= next16_name
            else:
                pass
        if next16_value in ['superkingdom','']:
            pass
        else:
            next17_= data[next16_==data['taxid']]['parent_taxid'].values[0]
            next17_name= names_sci[next17_==names_sci['taxid']]['name'].values[0]
            next17_value= data[next17_==data['taxid']]['rank'].values[0]
            if next17_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next17_value][i]= next17_name
            else:
                pass
        if next17_value in ['superkingdom','']:
            pass
        else:
            next18_= data[next17_==data['taxid']]['parent_taxid'].values[0]
            next18_name= names_sci[next18_==names_sci['taxid']]['name'].values[0]
            next18_value= data[next18_==data['taxid']]['rank'].values[0]
            if next18_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next18_value][i]= next18_name
            else:
                pass
        if next18_value in ['superkingdom','']:
            pass
        else:
            next19_= data[next18_==data['taxid']]['parent_taxid'].values[0]
            next19_name= names_sci[next19_==names_sci['taxid']]['name'].values[0]
            next19_value= data[next19_==data['taxid']]['rank'].values[0]
            if next19_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next19_value][i]= next19_name
            else:
                pass
        if next19_value in ['superkingdom','']:
            pass
        else:
            next20_= data[next19_==data['taxid']]['parent_taxid'].values[0]
            next20_name= names_sci[next20_==names_sci['taxid']]['name'].values[0]
            next20_value= data[next20_==data['taxid']]['rank'].values[0]
            if next20_value in ['species','genus','family','order','class','phylum','kingdom','superkingdom','clade']:
                blast_data[next20_value][i]= next20_name
            else:
                pass
        if next20_value in ['superkingdom','']:
            pass
        else:
            blast_data[next20_value][i]= 'Taxonomy_Sucks'
    else:
        blast_data['species'][i]= 'NA'
        blast_data['genus'][i]= 'NA'
        blast_data['family'][i]= 'NA'
        blast_data['order'][i]= 'NA'
        blast_data['class'][i]= 'NA'
        blast_data['phylum'][i]= 'NA'
        blast_data['kingdom'][i]= 'NA'
        blast_data['superkingdom'][i]= 'NA'
        blast_data['clade'][i]= 'NA'


blast_data.to_csv(output_csv,sep='\t',index=False)

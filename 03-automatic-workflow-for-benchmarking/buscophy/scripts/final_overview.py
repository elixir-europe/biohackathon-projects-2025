from os import listdir
from os.path import isfile, join
import os
import re
import subprocess
from Bio import Phylo
import csv
import sys 

def get_fasta_filenames(directory):

    onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]

    onlyfasta = [f for f in onlyfiles if re.search('.fasta', f)]

    return onlyfasta


def extract_taxa(directory):

    pattern='>'
    names=set()

    for file in get_fasta_filenames(directory):

        f = open(os.path.join(directory, file), "r")
        for line in f:
            if re.search(pattern, line):
                names.add(re.sub(pattern , '' , line.strip()))
                    #
    return sorted(names)

def build_header(directory):

    columns = ['File_name']
    
    for SPEC in extract_taxa(directory):
        columns.append(SPEC)
	
    columns.append('total_species')
    
    return columns 


def check_for_taxa(file, path,taxa):
    tax_dict=dict()
    for taxon in taxa:
        tax_dict[taxon]=0
    
    f = open(os.path.join(path, file), "r")
    
    for line in f:
        for taxon in taxa:
            if re.search(taxon, line):
                tax_dict[taxon]=1
    f.close()

    tax_list=list()
    for taxon in taxa:
        tax_list.append(tax_dict[taxon])

    total = sum(tax_list)
    tax_list.append(total)
    return tax_list

def main(IN_DIR,OUT_FILE ):
    header = build_header(IN_DIR)
    taxa=extract_taxa(IN_DIR)

    with open(OUT_FILE, 'w', buffering=1) as out_file:
        writer = csv.writer(out_file, delimiter='\t')

        writer.writerow(header)

        for file in get_fasta_filenames(IN_DIR):
            collection=[file]

            collection.extend(check_for_taxa(file, IN_DIR, taxa))

            writer.writerow(collection)


if __name__ == "__main__":
    #IN_DIR=snakemake.params.in_dir 
    #OUT_FILE=snakemake.output.table #output needs to be named in rule!
    IN_DIR=sys.argv[1]
    OUT_FILE=sys.argv[2]

    print(IN_DIR)
    print(type(OUT_FILE))

    main(IN_DIR,OUT_FILE )
    

            

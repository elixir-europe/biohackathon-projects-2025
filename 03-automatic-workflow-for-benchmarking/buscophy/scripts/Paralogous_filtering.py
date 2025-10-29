#!/usr/bin/env python3
import sys
from ete4 import PhyloTree
from ete4.smartview import Layout
import re
import os
from Bio import SeqIO
from Bio import AlignIO
import shutil
from Bio.Align import MultipleSeqAlignment

'''
(c) Code written by Nikita Kulikov(1), Iker Irisarri(2), Miguel Aranda(3); 
1) Leibniz Institute for the Analysis of Biodiversity Change (LIB)
2) Museo Nacional de Ciencias Naturales-CSIC (MNCN-CSIC)

This code loop through gene trees and determine speciation and duplication events. Then it removes inparalogous and outparalogous from the 
initial alignment, saving removed sequences in a separated file. It saves a phylogenetic tree in an extended newick format with the info
about speciation/duplication events. 


usage: python3 Paralogous_filtering.py  directory_with_trees directory_with_-corresponding_alignments

2025
'''

# First argument for a tree directory, second argument for a alignment directory
tree_file = sys.argv[1]
alignment_file =  sys.argv[2]

def extract_int_node_names(node):
    '''
    Rename nodes keeping just species names
    :param node: each node of the tree
    :return: renamed node
    '''
    name_parts = node.name.split('|')[0].split('_')
    if bool(re.search(r'\d+at\d+', name_parts[2])):
        my_node = '_'.join(name_parts[:2])
    else:
        my_node = '_'.join(name_parts[:3])

    return my_node

def load_and_root_tree(tfile, afile=None):
    '''
    Load and root a tree
    :param tfile: a tree file in a newick format
    :return: loaded and rooted tree
    '''
    if afile:
        t = PhyloTree(tfile, sp_naming_function=extract_int_node_names)
        t.link_to_alignment(afile)
    else:
        t = PhyloTree(tfile, sp_naming_function=extract_int_node_names)

    # Root the gene tree to midpoint
    out = t.get_midpoint_outgroup()
    if out is not t:
        t.set_outgroup(out)

    t.standardize()
    return t

def deal_with_inparalogs(my_inpar,al_f):
    '''
    From a provided list of inparalog sequences keeps the longest sequence only. It also creates a file with deleted
    inparaloug sequences.
    :param my_inpar, al_f: list of all inparalogousss sequences, alignment file
    :return: nothing
    '''
    records = list(SeqIO.parse(al_f, "fasta"))

    target_records = [rec for rec in records if rec.id in my_inpar]
    other_records = [rec for rec in records if rec.id not in my_inpar]
    removed_records = []
    removed_file = al_f.replace(".clipkit", ".clipkit.removed")

    if target_records:
        longest_target = max(target_records, key=lambda r: len(r.seq.replace("-","").replace("X","")))
        removed_records = [rec for rec in target_records if rec.id != longest_target.id]

        final_records = other_records + [longest_target]

        SeqIO.write(final_records, al_f, "fasta")

    if removed_records:
        with open(removed_file, "a") as rem_f:
            SeqIO.write(removed_records, rem_f, "fasta")


def count_number_of_sp(child):
    '''
    This function extracts set of species from the outparalogous branches
    :param child: one of the branches
    :return: set of species of this branch
    '''
    processed_nodes = set()
    for node_name in child:
        name_parts = node_name.split('|')[0].split('_')
        if len(name_parts) > 2 and re.search(r'\d+at\d+', name_parts[2]):
            my_node = '_'.join(name_parts[:2])
        else:
            my_node = '_'.join(name_parts[:3])
        processed_nodes.add(my_node)

    return processed_nodes

def keep_on_branch(child1,child2):
    '''
    This function extracts set of species from the outparalogous branches
    :param child1, child2: two branches with outparalogous
    :return: which branch has more species
    '''
    num_of_sp1 = count_number_of_sp(child1)
    num_of_sp2 = count_number_of_sp(child2)

    if len(num_of_sp1) > len(num_of_sp2):
        return "child1"
    else:
        return "child2"

def deal_with_outparalogs(child1,child2,alignment_filename):
    '''
    From provided branches with outparaloug sequences keeps the branch with the highest number of sequences. It also creates a file with deleted
    outparaloug sequences.
    :param child1, child2, alignment_filename: first child branch of  a duplication event, second branch, alignment file
    :return: nothing
    '''

    which_br_to_keep = keep_on_branch(child1,child2)
    alignment = AlignIO.read(alignment_filename, "fasta")

    removed_records = []
    removed_file = alignment_filename.replace(".clipkit", ".clipkit.removed")
    if which_br_to_keep == "child1":
        filtered_records = [record for record in alignment if record.id not in child2]
        removed_records = [record for record in alignment if record.id in child2]
    else:
        filtered_records = [record for record in alignment if record.id not in child1]
        removed_records = [record for record in alignment if record.id in child1]

    filtered_alignment = MultipleSeqAlignment(filtered_records)

    AlignIO.write(filtered_alignment, alignment_filename, "fasta")

    if removed_records:
        with open(removed_file, "a") as rem_f:
            SeqIO.write(removed_records, rem_f, "fasta")


def detect_duplication_events(t,al_f):
    '''
    Modifies a tree by determining duplication and speciation events . It also call for functions to remove in- and outparalogous from alignments
    :param t, al_f: tree file, alignment file
    :return: nothing
    '''
    events = t.get_descendant_evol_events()


    for ev in events:

        r = {'S': 'Orthology', 'D': 'Paralogy'}[ev.etype]  # relationship

        if ev.etype == "D":
            print("We found a paralogy!")

            first_br = count_number_of_sp(ev.in_seqs)
            second_br = count_number_of_sp(ev.out_seqs)

            if len(first_br) == 1 and first_br == second_br:
                print("Inparalogy was detected!")
                my_inparalogy = ev.in_seqs.union(ev.out_seqs)
                deal_with_inparalogs(my_inparalogy,al_f)

            else:
                print("Outparalogy was detected!")
                deal_with_outparalogs(ev.in_seqs,ev.out_seqs,al_f)

        else:

            print("Orthology was detected!")

def view_tree(tfile, afile):

    # get the lists of all files in the directories
    tree_list = os.listdir(tfile)
    alignment_list = os.listdir(afile)

    for tree_file in tree_list:
        # look for the treefile extension
        if tree_file.endswith(".treefile"):
            for alignment_file in alignment_list:
                # check if the alignment file corresponds to the tree file
                tree_file_check = tree_file.split(".")[0].split("_")
                tree_file_check_str =  '_'.join(tree_file_check[:2])
                alignment_file_check = alignment_file.split(".")[0].split("_")
                alignment_file_check_str =  '_'.join(alignment_file_check[:2])
                if alignment_file_check_str == tree_file_check_str:
                    tree_file = os.path.join(tfile, tree_file)
                    alignment_file = os.path.join(afile, alignment_file)
                    # load and root trees to a midpoint using ete4
                    t = load_and_root_tree(tree_file, alignment_file)
                    # preserve the original alignment file
                    alignment_file_copy = alignment_file.replace(".clipkit", ".clipkit.original")
                    shutil.copy(alignment_file, alignment_file_copy)
                    # detecte speciation/duplication events and remove in- and outparalogous from the alignment
                    detect_duplication_events(t,alignment_file)
                    # save extended newick trees with the duplication/speciation info
                    modified_tree_filename = tree_file.replace(".treefile",".extended_treefile")
                    t.write(modified_tree_filename,props=['evoltype'])


#excute only if called from the command line
if __name__ == '__main__':
        view_tree(sys.argv[1], sys.argv[2])


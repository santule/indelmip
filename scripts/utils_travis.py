import pickle
import numpy as np
from itertools import accumulate
import operator
from pysam import FastaFile,FastxFile
from os import walk
from ete3 import Tree
import os
import sys
import glob
import numpy as np
import pandas as pd


def total_and_length(fasta_file):
    all_seq = []
    total_sequences = 0
    with FastxFile(fasta_file) as fh:
        for entry in fh:
            # add start and end string to the sequence
            seq_name = entry.sequence
            seq_binary   = [0 if s == '-' else 1 for s in seq_name]
            all_seq.append(seq_binary)
            total_sequences += 1
            sequence_length = len(seq_binary)
            
    return all_seq,sequence_length,total_sequences
                
def eff_pos_stat(all_seq,total_sequences,sequence_length):
    
    sum_all_pos = (np.array(all_seq).sum(axis=0) == total_sequences) + 0
    # count consecutive 1s of size > 2
    ones_sets  = []
    cons_ones  = []
    found_one  = 0

    for ind,i in enumerate(sum_all_pos):
        if i == 1:
            if found_one == 1:
                cons_ones.append(ind)
            else:
                cons_ones.append(ind)
                found_one = 1
        else:
            found_one = 0
            if len(cons_ones) >= 2:
                ones_sets.append(cons_ones)
            cons_ones = []

    # calculate effective positions
    re_pos = 0
    for o in ones_sets:
        re_pos += len(o) - 1

    eff_pos = sequence_length - re_pos
    eff_pos_percent = round(eff_pos/(eff_pos + re_pos) * 100,2)

    return eff_pos_percent

def avg_edges_per_pos(ex_pickle_file,ancestor_pickle_file):
    with open(ancestor_pickle_file, 'rb') as f:
        ancestor_data = pickle.load(f)

    with open(ex_pickle_file, 'rb') as f:
        extant_data = pickle.load(f)
        
    sequence_length = len(list(extant_data.values())[0])
    edges_out = ancestor_data[1]
    edges_in  = ancestor_data[2]
    
    eo_len_list = [len(eo) for eo in list(edges_out.values())]
    ei_len_list = [len(ei) for ei in list(edges_in.values())]
    
    
    all_edges = list( map(operator.add, eo_len_list, ei_len_list) )     
    avg_edges  = sum(all_edges) / len(all_edges)

    return avg_edges

def get_stats(sub_folder):
    
    ex_file   = '/extants.aln'
    tree_file = '/input_tree.nwk'
    all_seq_file  = '/all_sequences.aln'
    ex_fasta_file = sub_folder + ex_file
    ex_pickle_file = sub_folder + '/extant_data.pkl'
    ancestor_pickle_file = sub_folder + '/ancestor_info.pkl'
    
    all_seq,sequence_length,total_sequences = total_and_length(ex_fasta_file)
    eff_pos_percent = eff_pos_stat(all_seq,total_sequences,sequence_length)
    avg_edges = avg_edges_per_pos(ex_pickle_file,ancestor_pickle_file)

    return sequence_length,total_sequences,eff_pos_percent,avg_edges

# seperate ancestor from extants
def seperate_ancestor_extants(folder,trim_length):
    tree_file = open(folder + 'input_tree.nwk',"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)
    sequences_fasta_info = FastaFile(folder + 'all_sequences.aln')

    with open(folder + 'extants.aln',mode='w') as fout:
        for n in tree.traverse(): # level order
            if n.is_leaf():
                extant_sequence = sequences_fasta_info.fetch(n.name)
                extant_name = n.name
                fout.write('>' + str(extant_name) + '\n')
                if trim_length != 0:
                    extant_sequence = extant_sequence[0:trim_length]
                fout.write(str(extant_sequence) + '\n')

    # remove *.fai file
    os.remove(folder + 'all_sequences.aln.fai')

# stats on travis data
def check_stats(folder,detail_print=False):
    # get all files in the folder
    all_ex_aln_files = []
    ex_file   = 'extants.aln'
    tree_file = 'input_tree.nwk'
    all_seq_file  = 'all_sequences.aln'

    for (sub_folder, _, _) in walk(folder):
        if sub_folder != folder:
            sub_folder = sub_folder + '/'
            print(f"Processing synthetic indel tree {sub_folder}")

            # seperate extants from ancestors
            seperate_ancestor_extants(sub_folder,tree_file,all_seq_file)

            # check stats
            all_seq = []
            total_sequences = 0
            with FastxFile(sub_folder + ex_file) as fh:
                for entry in fh:
                    # add start and end string to the sequence
                    seq_name = entry.sequence
                    seq_binary   = [0 if s == '-' else 1 for s in seq_name]
                    all_seq.append(seq_binary)
                    total_sequences += 1
                    sequence_length = len(seq_binary)

            print(f"Sequence length {sequence_length}")
            print(f"Total Sequences {total_sequences}")

            sum_all_pos = (np.array(all_seq).sum(axis=0) == total_sequences) + 0

            # count consecutive 1s of size > 2
            ones_sets  = []
            cons_ones  = []
            found_one  = 0

            for ind,i in enumerate(sum_all_pos):
                if i == 1:
                    if found_one == 1:
                        cons_ones.append(ind)
                    else:
                        cons_ones.append(ind)
                        found_one = 1
                else:
                    found_one = 0
                    if len(cons_ones) >= 2:
                        ones_sets.append(cons_ones)
                    cons_ones = []

            # calculate effective positions
            re_pos = 0
            for o in ones_sets:
                re_pos += len(o) - 1

            eff_pos = sequence_length - re_pos
            eff_pos_percent = round(eff_pos/(eff_pos + re_pos) * 100,2)

            print(f"redundant positions:{re_pos}")
            eff_pos = sequence_length - re_pos
            print(f"effective positions:{eff_pos}")
            print(f"effective positions %:{round(eff_pos/(eff_pos + re_pos) * 100,2)}")


if __name__ == "__main__":
    folder = sys.argv[1] 
    trim_length = int(sys.argv[2])
    #check_stats(folder,detail_print=False)
    seperate_ancestor_extants(folder,trim_length)

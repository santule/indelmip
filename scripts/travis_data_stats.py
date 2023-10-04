import pickle
import numpy as np
from itertools import accumulate
import operator
from pysam import FastaFile,FastxFile
from os import walk
from ete3 import Tree
import os
import sys

# seperate ancestor from extants
def seperate_ancestor_extants(folder,tree_file,all_seq_file):
    tree_file = open(folder + tree_file,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)
    sequences_fasta_info = FastaFile(folder + all_seq_file)

    with open(folder + 'extants.aln',mode='w') as fout:
        for n in tree.traverse(): # level order
            if n.is_leaf():
                extant_sequence = sequences_fasta_info.fetch(n.name)
                extant_name = n.name
                fout.write('>' + str(extant_name) + '\n')
                fout.write(str(extant_sequence) + '\n')

    # remove *.fai file
    os.remove(folder + 'all_sequences.aln.fai')

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
  folder = sys.argv[1] #'/Users/sanjanatule/Documents/uq/Projects/Indels/indelmip/data/travis_500/'
  check_stats(folder,detail_print=False)

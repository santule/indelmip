import pickle
import numpy as np
from itertools import accumulate
import operator
from pysam import FastaFile,FastxFile
from os import walk
from ete3 import Tree
import os

def check_stats(folder,detail_print=False):
    # get all files in the folder
    all_ex_aln_files = []
    ex_file   = 'extants.aln'
    tree_file = 'input_tree.nwk'
    all_seq_file  = 'all_sequences.aln'

    for (sub_folder, _, _) in walk(folder):
        if sub_folder != folder:
            sub_folder = sub_folder + '/'
            print(sub_folder)
            desired_length  = int(str(sub_folder).split('/')[-2].split('l')[-1].split('e')[0])
            desired_eff     = int(str(sub_folder).split('/')[-2].split('l')[-1].split('e')[1])
            desired_num_seq = int(str(sub_folder).split('/')[-2].split('l')[0].split('t')[1])

            print(f"Processing synthetic indel tree {sub_folder}")

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

                    if sequence_length != desired_length:
                        print("Sequence length not correct.Failed.")

            if total_sequences != desired_num_seq:
                print("Total extant sequences are incorrect.Failed.")

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

            if eff_pos_percent <= desired_eff + 5 and eff_pos_percent >= desired_eff - 5:
                print("Successfully created")
                if detail_print:
                    print(f"********** STATS **********")
                    print(f"file name:{sub_folder + ex_file}")
                    print(f"total_sequences:{total_sequences}")
                    print(f"sequence_length:{sequence_length}")
                    print(f"redundant positions:{re_pos}")
                    eff_pos = sequence_length - re_pos
                    print(f"effective positions:{eff_pos}")
                    print(f"effective positions %:{round(eff_pos/(eff_pos + re_pos) * 100,2)}")
            else:
                print("Eff positions not correct.Failed")
                print(f"********** STATS **********")
                print(f"file name:{sub_folder + ex_file}")
                print(f"total_sequences:{total_sequences}")
                print(f"sequence_length:{sequence_length}")
                print(f"redundant positions:{re_pos}")
                eff_pos = sequence_length - re_pos
                print(f"effective positions:{eff_pos}")
                print(f"effective positions %:{round(eff_pos/(eff_pos + re_pos) * 100,2)}")

if __name__ == "__main__":
  folder = '/Users/sanjanatule/Documents/uq/Projects/Indels/indelmip/data/travis/'
  check_stats(folder,detail_print=False)

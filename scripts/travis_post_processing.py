import pickle
import numpy as np
from itertools import accumulate
import operator
from pysam import FastaFile,FastxFile
from os import walk
from ete3 import Tree
import os

# rectify effective pos
def check_rectify_eff(fasta_file,desired_eff_pos):

    all_seq = []
    total_sequences = 0
    with FastxFile(fasta_file) as fh:
        for entry in fh:
            # add start and end string to the sequence
            seq_name = entry.sequence
            seq_binary   = [0 if s == '-' else 1 for s in seq_name]
            all_seq.append(seq_binary)
            total_sequences += 1
            sequence_length = len(seq_name)

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


    if eff_pos_percent <= desired_eff_pos + 5 and eff_pos_percent >= desired_eff_pos - 5:

        print(f"********** FINISHED **********")
        print("Created synthetic indels successfully!!!")
        print(f"file name:{fasta_file}")
        print(f"total_sequences:{total_sequences}")
        print(f"sequence_length:{sequence_length}")
        print(f"redundant positions:{re_pos}")
        eff_pos = sequence_length - re_pos
        print(f"effective positions:{eff_pos}")
        print(f"effective positions %:{round(eff_pos/(eff_pos + re_pos) * 100,2)}")
        print(f"*****************************")

        return 1
    return 0

# seperate ancestor from extants
def seperate_ancestor_extants(folder,tree_file,all_seq_file,desired_length,start_pos):
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
                fout.write(str(extant_sequence[start_pos:start_pos + desired_length]) + '\n')
    # remove *.fai file
    os.remove(folder + 'all_sequences.aln.fai')

# make extants
def make_extants(folder):
    
    # get all files in the folder
    all_ex_aln_files = []
    ex_file       = 'extants.aln'
    tree_file     = 'input_tree.nwk'
    all_seq_file  = 'all_sequences.aln'

    for (sub_folder, _, _) in walk(folder):
        if sub_folder != folder:
            sub_folder = sub_folder + '/'
            desired_length = int(str(sub_folder).split('/')[-2].split('l')[-1].split('e')[0])
            desired_eff    = int(str(sub_folder).split('/')[-2].split('l')[-1].split('e')[1])
            print(f"Tree: {sub_folder}")

            for start_pos in range(0,3554 - desired_length):
                seperate_ancestor_extants(sub_folder,tree_file,all_seq_file,desired_length,start_pos) # seperate ancestors from extants
                if check_rectify_eff(sub_folder + ex_file,desired_eff):
                    break

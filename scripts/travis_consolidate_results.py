''' python script to consolidate the results from travis experiment '''
import glob, os
from os import walk
from pysam import FastaFile,FastxFile
from ete3 import Tree
import numpy as np
import pickle
import operator
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


def csv_pivot(csv_file,value_col_name):
    pd_data = pd.read_csv(csv_file)
    print(pd_data.head(5))
    pd_data = pd_data.pivot(index=['protein_family','sequence_length','total_sequences','eff_pos_percent','avg_edges'], columns=['method'],values=value_col_name).reset_index()
    # pd_data.columns = pd_data.columns.droplevel()
    # pd_data = pd_data.reset_index()
    pd_data.to_csv(csv_file)
    print("finished summarizing experimental results")

### MAIN PROGRAM ###

synthetic_stat = {}
sample_file = 'extants.aln'
experiment_result_folder = '/media/WorkingSpace/Share/mipindel/evaluation/synthetic/'
data_folder = '/media/WorkingSpace/Share/mipindel/data/travis'

# get the stat of synthetic dataset
for (sub_folder, _, _) in walk(data_folder):
    os.chdir(sub_folder)
    if os.path.exists(sample_file):
        print(sub_folder)
        sequence_length,total_sequences,eff_pos_percent,avg_edges = get_stats(sub_folder) 
        data_idx = sub_folder.split('/')[-1]
        synthetic_stat[data_idx] = (sequence_length,total_sequences,round(eff_pos_percent,2),round(avg_edges,2))


# collect experiment results
print('Getting Indel scores')
with open(experiment_result_folder + 'indel_diff_evaluation.csv', 'w') as outfile:
    outfile.write('protein_family,sequence_length,total_sequences,eff_pos_percent,avg_edges,method,score\n')
    for (sub_folder, _, _) in walk(data_folder):
        os.chdir(sub_folder)
        for fname in glob.glob("*indscore*"):
            with open(fname) as infile:
                for line in infile:
                    pr = sub_folder.split('/')[-1]
                    outfile.write(pr + ',' + str(synthetic_stat[pr][0]) + ',' + str(synthetic_stat[pr][1]) + ',' + str(synthetic_stat[pr][2]) + ',' + str(synthetic_stat[pr][3]) + ',' + line + '\n' )

csv_pivot(experiment_result_folder + 'indel_diff_evaluation.csv','score')

print('Getting Out of Distribution scores')
with open(experiment_result_folder + 'out_of_dist_evaluation.csv', 'w') as outfile:
    outfile.write('protein_family,sequence_length,total_sequences,eff_pos_percent,avg_edges,method,out_of_dist_percent\n')
    for (sub_folder, _, _) in walk(data_folder):
        os.chdir(sub_folder)
        for fname in glob.glob("*out_dist_percent*"):
            with open(fname) as infile:
                for line in infile:
                    pr = sub_folder.split('/')[-1]
                    outfile.write(pr + ',' + str(synthetic_stat[pr][0]) + ',' + str(synthetic_stat[pr][1]) + ',' + str(synthetic_stat[pr][2]) + ',' + str(synthetic_stat[pr][3]) + ',' + line + '\n' )
csv_pivot(experiment_result_folder + 'out_of_dist_evaluation.csv','out_of_dist_percent')

print('Getting 3 Mutation scores')
with open(experiment_result_folder + 'ancestors_with_3_mutation.csv', 'w') as outfile:
    outfile.write('protein_family,sequence_length,total_sequences,eff_pos_percent,avg_edges,method,ancestors_3_mutation_percent\n')
    for (sub_folder, _, _) in walk(data_folder):
        os.chdir(sub_folder)
        for fname in glob.glob("*ancestors_with_3_mut*"):
            with open(fname) as infile:
                for line in infile:
                    pr = sub_folder.split('/')[-1]
                    outfile.write(pr + ',' + str(synthetic_stat[pr][0]) + ',' + str(synthetic_stat[pr][1]) + ',' + str(synthetic_stat[pr][2]) + ',' + str(synthetic_stat[pr][3]) + ',' + line + '\n' )
csv_pivot(experiment_result_folder + 'ancestors_with_3_mutation.csv','ancestors_3_mutation_percent')

print('Getting 2 Mutation scores')
with open(experiment_result_folder + 'ancestors_2_mutation.csv', 'w') as outfile:
    outfile.write('protein_family,sequence_length,total_sequences,eff_pos_percent,avg_edges,method,ancestors_2_mutation_percent\n')
    for (sub_folder, _, _) in walk(data_folder):
        os.chdir(sub_folder)
        for fname in glob.glob("*ancestors_with_2_mut*"):
            with open(fname) as infile:
                for line in infile:
                    pr = sub_folder.split('/')[-1]
                    outfile.write(pr + ',' + str(synthetic_stat[pr][0]) + ',' + str(synthetic_stat[pr][1]) + ',' + str(synthetic_stat[pr][2]) + ',' + str(synthetic_stat[pr][3]) + ',' + line + '\n' )
csv_pivot(experiment_result_folder + 'ancestors_2_mutation.csv','ancestors_2_mutation_percent')

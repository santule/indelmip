from pysam import FastaFile,FastxFile
from ete3 import Tree
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import getopt

''' function to get the sequence levels in the tree '''
def get_sequence_level_in_tree(nwk_file_path):
    tree_file = open(nwk_file_path,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)

    sequence_tree_level = {}
    ancestor_names = []
    extant_names = []

    # get all ancestor and extant names
    for n in tree.traverse():
        if not n.is_leaf():
            ancestor_names.append(n.name)
        else:
            extant_names.append(n.name)

    # traverse the tree and add levels.
    level = 0
    for n in tree.traverse():
        if n.up is not None: # root node
            n.add_features(level = n.up.level + 1)
            sequence_tree_level[n.name] = n.up.level + 1
        else:
            n.add_features(level = level)
            sequence_tree_level[n.name] = level
    return sequence_tree_level,extant_names,ancestor_names

''' function to build common pattern vocab '''
def create_pattern_data(indel_file,method_name,extant_names,nwk_file_path):
    # get the indel patterns and give them unique number
    sequence_pattern_dict  = {} # common dictionary of all indel patterns shared by  a solution
    pattern_type = {} # pattern id - extant or ancestor
    sequence_name_pattern_idx = {} # seq name - pattern id
    extant_patterns = [] # get all extant pattern
    pattern = 1
    sequences_fasta_info = FastaFile(indel_file)
    tree_file = open(nwk_file_path,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)

    for n in tree.traverse():
        seq_name = n.name
        indel_sequence = sequences_fasta_info.fetch(n.name)

        if method_name == 'mip':
            indel_sequence = indel_sequence[1:-1] # remove start and end node

        if seq_name in extant_names:
            extant_patterns.append(indel_sequence)
        else:
            if sequence_pattern_dict.get(indel_sequence):
                sequence_name_pattern_idx[seq_name]   = sequence_pattern_dict[indel_sequence]
            else:
                sequence_pattern_dict[indel_sequence] = pattern
                sequence_name_pattern_idx[seq_name]   = pattern
                pattern_type[pattern] = 'A' # default allocation, to be changed later
                pattern = pattern + 1

    #### Classify patterns into new / old ####
    for k,v in sequence_pattern_dict.items():
        if k in extant_patterns: # if the pattern is in extant pattern
            pattern_type[v] = 'E'
    return pattern_type,sequence_name_pattern_idx

''' function to combine data in dataframe'''
def combine_data_df(sequence_tree_level,sequence_name_pattern_idx,pattern_type):

    # wrangle data for visualisation
    indel_visual_df1 = pd.DataFrame(sequence_tree_level.items(),columns = ['name','tree_level'])
    indel_visual_df2 = pd.DataFrame(sequence_name_pattern_idx.items(),columns = ['name','pattern'])
    indel_visual_df3 = pd.DataFrame(pattern_type.items(),columns = ['pattern','type'])
    indel_visual_df  = pd.merge(indel_visual_df2, indel_visual_df1, how= "inner" , on="name")
    indel_visual_df  = pd.merge(indel_visual_df, indel_visual_df3,  how= "inner" , on="pattern")
    indel_visual_df  = indel_visual_df.drop(columns='name')
    indel_visual_df['dummy_value'] = np.where(indel_visual_df['type'] =='E', 1,2)
    indel_visual_df = indel_visual_df.drop(['type'],axis=1)
    indel_visual_df  = indel_visual_df.drop_duplicates()
    indel_visual_df  = indel_visual_df.pivot(index='pattern', columns= 'tree_level',values = 'dummy_value')
    indel_visual_df  = indel_visual_df.fillna(0)
    indel_visual_df  = indel_visual_df.astype(int)
    return indel_visual_df

''' function to create heatmap'''
def create_plot(indel_visual_df,output_file_name,method):
    # heatmap
    plt.rcParams["figure.figsize"] = [10,10]
    colors = ["#F5F3EE", "#738CB1","#77C24D"] # 1 - extant (blue) , 2 - ancestor(green)
    sns_pp = sns.heatmap(indel_visual_df,cmap = colors,cbar=False)
    sns_pp.set(xlabel='Tree Level' , ylabel='Indel Pattern')
    plt.title(f"Visualisation for {method}")
    # save
    plt.savefig(output_file_name)

# main function
def visualise_solution(nwk_file_path,indel_file,method):
    # tree levels
    sequence_tree_level,extant_names,ancestor_names = get_sequence_level_in_tree(nwk_file_path)
    # patterns
    pattern_type,sequence_name_pattern_idx = \
    create_pattern_data(indel_file,method,extant_names,nwk_file_path)
    # create df
    visual_df = combine_data_df(sequence_tree_level,sequence_name_pattern_idx,pattern_type)
    # plots
    create_plot(visual_df,method + '_pattern.jpeg',method)
    print(f"Output file created :{method}_pattern.jpeg")

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python visualise_sol.py -i indel fasta file | -n newick tree file')
    exit()

# main function
def main(indel_fasta_solution_file,nwk_file_path,method):
    # get the parimony score
    visualise_solution(indel_fasta_solution_file,nwk_file_path,method)

if __name__ == "__main__":
  main(nwk_file_path,indel_fasta_solution_file,method)

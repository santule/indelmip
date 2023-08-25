''' script to score indel output and calculate the different indel events '''

from ete3 import Tree
from pysam import FastaFile,FastxFile
import sys
import getopt

''' Function to calculate parsimony score for the whole tree for a given indel solution'''
def sequence_distance_score(str1,str2):
    dis = 0
    prev_dis = 0

    for i in range(0,len(str1)):
        if str1[i] != str2[i]:  # not matching
            if prev_dis == 0:   # previous matches
                dis += 6
                prev_dis = 1
            else:
                dis += 1
        else:
            prev_dis = 0
    return dis

def score_tree_indels(treefile,indelfastafile):
    parsiscore = 0
    # load the fasta file
    indel_pattern = FastaFile(indelfastafile)
    # load the tree
    tree_file = open(treefile,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)

    # load pattern on the tree
    for n in tree.traverse():
        if n.is_leaf() == False:
            current_node = n.name
            current_node_sequence = indel_pattern.fetch(current_node)
            child_seq_1 = indel_pattern.fetch(n.children[0].name)
            child_seq_2 = indel_pattern.fetch(n.children[1].name)

            # calculate score
            parsiscore += sequence_distance_score(current_node_sequence,child_seq_1)
            parsiscore += sequence_distance_score(current_node_sequence,child_seq_2)

    return parsiscore

# main function
def main(nwk_file_path,indel_fasta_solution_file):
    print(indel_fasta_solution_file)
    # get the parimony score
    parsiscore = score_tree_indels(nwk_file_path,indel_fasta_solution_file)
    print(f"The overall indel score is {parsiscore}")

if __name__ == "__main__":
  main(nwk_file_path,indel_fasta_solution_file)

''' script to calculate objective score for a solution '''

from ete3 import Tree
from pysam import FastaFile,FastxFile
import sys
import getopt

''' Function to calculate objective score for the whole tree for a given indel solution'''
def sequence_obj_score(str1,str2):
    objscore = 0
    prev_dis = 0
    curr_dis = 0

    for i in range(0,len(str1)):
        curr_dis  = int(str1[i]) - int(str2[i])
        
        if curr_dis != 0 and curr_dis != prev_dis:
            objscore += 2
            
        objscore += abs(curr_dis)
        prev_dis = curr_dis
        
    return objscore

def score_tree_objective(treefile,indelfastafile):
    obj_score = 0
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
            obj_score += sequence_obj_score(current_node_sequence,child_seq_1)
            obj_score += sequence_obj_score(current_node_sequence,child_seq_2)

    return obj_score

# main function
def main(nwk_file_path,indel_fasta_solution_file):
    # get the parimony score
    obj_score = score_tree_objective(nwk_file_path,indel_fasta_solution_file)
    print(f"The objective score is {obj_score}")
    return obj_score

if __name__ == "__main__":
    obj_score = main(nwk_file_path,indel_fasta_solution_file)

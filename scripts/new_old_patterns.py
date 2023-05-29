''' script to score indel output and calculate the different indel events '''

from ete3 import Tree
from pysam import FastaFile,FastxFile
import sys
import getopt


def score_tree_indels(treefile,indelfastafile):
    parsiscore = 0
    # load the fasta file
    indel_sequences = FastaFile(indelfastafile)
    # load the tree
    tree_file = open(treefile,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)
    existing_extant_patterns = []

    # load extant pattern on the tree
    for n in tree.traverse():
        if n.is_leaf() == True:
            indel = indel_sequences.fetch(n.name)
            existing_extant_patterns.append(indel)

    # count old and new patterns
    extant_pattern = 0
    new_pattern = 0
    for n in tree.traverse():
        if n.is_leaf() == False:
            if indel_sequences.fetch(n.name) in existing_extant_patterns:
                extant_pattern += 1
            else:
                new_pattern +=1

    return extant_pattern,new_pattern


# main function
def main(nwk_file_path,indel_fasta_solution_file):
    # get the parimony score
    extant_pattern,new_pattern = score_tree_indels(nwk_file_path,indel_fasta_solution_file)
    print(f"Total {extant_pattern} extant patterns and {new_pattern} new patterns used to create ancestors")

if __name__ == "__main__":
  main(nwk_file_path,indel_fasta_solution_file)

''' script to score indel output and calculate the different indel events '''

from ete3 import Tree
from pysam import FastaFile,FastxFile
import sys
import getopt
from operator import add

''' function to calculate the total indel events at each node of the tree'''
def count_mutations(p,n,c1,c2):
    i3 = 0
    for i in range(0,len(n)):
        if n[i] != p[i] and n[i] != c1[i] and n[i] != c2[i]:
            i3 += 1

    if i3 > 0:
        return 1
    else:
        return 0

def count_indel_events(treefile,indelfastafile):

    total_ancestors = 0
    total_i3_mut = 0

    # load the tree
    tree_file = open(treefile,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)

    # traverse tree
    indel_pattern = FastaFile(indelfastafile)

    for n in tree.traverse():
        if n.is_leaf() == False:
            total_ancestors += 1
            current_node = n.name
            current_node_sequence = indel_pattern.fetch(current_node)
            child_seq_1 = indel_pattern.fetch(n.children[0].name)
            child_seq_2 = indel_pattern.fetch(n.children[1].name)

            if n.is_root() == False: # if root then parent sequence same as the root
                parent_node  = n.up.name
                parent_node_sequence = indel_pattern.fetch(parent_node)
            else:
                parent_node_sequence = current_node_sequence

            i3_mutations  = count_mutations(parent_node_sequence,current_node_sequence,child_seq_1,child_seq_2)
            total_i3_mut +=  i3_mutations

    return total_i3_mut , round((total_i3_mut / total_ancestors)* 100,2)


# main function
def main(nwk_file_path,indel_fasta_solution_file):

    # get the parimony score
    total_i3_mut,percent_ancestor_with_3_mut = count_indel_events(nwk_file_path,indel_fasta_solution_file)
    print(total_i3_mut,percent_ancestor_with_3_mut)
    print(f"Total ancestors with pattern 3 mutations away:{total_i3_mut}")
    print(f"% ancestors with pattern 3 mutations away:{percent_ancestor_with_3_mut}")
    return total_i3_mut,percent_ancestor_with_3_mut

if __name__ == "__main__":
  total_i3_mut,percent_ancestor_with_3_mut = main(nwk_file_path,indel_fasta_solution_file)

''' script to score indel output and calculate the different indel events '''

from ete3 import Tree
from pysam import FastaFile,FastxFile
import sys
import getopt
from operator import add

''' function to calculate the total indel events at each node of the tree'''
def count_mutations(p,n,c1,c2):
    tn,i1a,i1b,i2a,i2b,i3 = 0,0,0,0,0,0

    for i in range(0,len(n)):
        if n[i] == p[i] == c1[i] == c2[i]:
            continue
        elif n[i] == p[i]: #(node and parent position are equal)
            if c1[i] == c2[i] and c1[i] != n[i]:
                i2b += 1
            if c1[i] != c2[i]: # kids are not equal
                i1b += 1
        elif n[i] != p[i]: #(node not equal to parent)
            if c1[i] != c2[i]:
                i2a += 1
            if c1[i] == c2[i] and c1[i] == n[i]:
                i1a += 1
            if c1[i] == c2[i] and c1[i] != n[i]:
                i3 += 1
    if i1a + i1b + i2a + i2b + i3 == 0: # no mutations
        tn = 0
    else:
        tn = 1
    return [tn,i1a,i1b,i2a,i2b,i3]

def count_indel_events(treefile,indelfastafile):
    total_mut = [0,0,0,0,0,0] #event flag,i1a,i1b,i2a,i2b,i3

    # load the tree
    tree_file = open(treefile,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)

    # traverse tree
    indel_pattern = FastaFile(indelfastafile)
    total_ancestors = 0
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
        mut_ret   = count_mutations(parent_node_sequence,current_node_sequence,child_seq_1,child_seq_2)
        total_mut =  list(map(add, total_mut, mut_ret))
    return total_mut, total_ancestors * len(parent_node_sequence)


# main function
def main(nwk_file_path,indel_fasta_solution_file):

    # get the parimony score
    indel_events_list,total_positions = count_indel_events(nwk_file_path,indel_fasta_solution_file)
    print(indel_events_list)
    print(f"Total sequences different from its neighbors in 1 or more positions:{indel_events_list[0]}")
    print(f"Total positions 1 mutation away from its parent:{indel_events_list[1]}")
    print(f"Total positions 1 mutation away from either of its children:{indel_events_list[2]}")
    print(f"Total positions 1 mutation away from parent and 1 mutation away from one of its children:{indel_events_list[3]}")
    print(f"Total positions 2 mutations away from its children:{indel_events_list[4]}")
    print(f"Total positions 3 mutations away:{indel_events_list[5]}")
    print(f"Total positions {total_positions}")

    return (indel_events_list[5])/(total_positions), indel_events_list[5]

if __name__ == "__main__":
  main(nwk_file_path,indel_fasta_solution_file)

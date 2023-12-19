''' script to score indel event difference for pairwise ancestor-child across the whole tree '''

from ete3 import Tree
from pysam import FastaFile,FastxFile

''' Function to calculate indel events for the whole tree for a given indel solution'''
def indel_events(str1,str2):
    dis = 0
    prev_dis = 0

    for i in range(0,len(str1)):
        curr_dis = int(str1[i]) - int(str2[i])
        
        if curr_dis != 0 and curr_dis != prev_dis:
            dis += 1
        prev_dis = curr_dis
        
    return dis

def score_tree_indels(treefile,indelfastafile):
    indelevents_cnt = 0
    
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
            indelevents_cnt += indel_events(current_node_sequence,child_seq_1)
            indelevents_cnt += indel_events(current_node_sequence,child_seq_2)

    return indelevents_cnt

# main function
def score_fn(nwk_file_path,indel_fasta_solution_file):
    # get the parimony score
    indelevents = score_tree_indels(nwk_file_path,indel_fasta_solution_file)
    print(f"The total indel events are {indelevents}")
    return indelevents

if __name__ == "__main__":
    indelevents = score_fn(nwk_file_path,indel_fasta_solution_file)

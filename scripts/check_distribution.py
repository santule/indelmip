''' script to check if ancestor indel pattern comes from extant indel distribution '''

from ete3 import Tree
from pysam import FastaFile,FastxFile
from collections import defaultdict

#### helper functions
def next_pos(str1,curr_pos,seq_len,gap_char = '-'):
    start_pos = curr_pos + 1

    while(start_pos < len(str1)):
        if str1[start_pos] != gap_char:
            return start_pos
        else:
            start_pos = start_pos + 1
    return seq_len

def convert_to_edges(seq_str,seq_pog_dict):
    seq_len = len(seq_str)
    ind = 0
    while(ind < seq_len - 1):
        if seq_str[ind] != '-':
            curr_ind = ind
            ind = next_pos(seq_str,curr_ind,seq_len - 1) # find the next filled position
            seq_pog_dict[curr_ind].add(ind)
        else:
            ind = ind + 1
    return seq_pog_dict

def create_extant_distribution(input_file):
    seq_pog_dict = defaultdict(set)
    with FastxFile(input_file) as fh:
        for entry in fh:
            # add start and end string to the sequence
            seq_name = entry.name
            new_sequence = 'x' + entry.sequence + 'x'
            sequence_length = len(new_sequence)

            # convert to edges
            seq_pog_dict = convert_to_edges(new_sequence,seq_pog_dict)

    return seq_pog_dict

def convert_indel_edge(seq_str):
    a_seq_pog_dict = {}
    seq_len = len(seq_str)
    ind = 0
    while(ind < seq_len - 1):
        if seq_str[ind] != '-':
            curr_ind = ind
            ind = next_pos(seq_str,curr_ind,seq_len - 1,'0') # find the next filled position
            a_seq_pog_dict[curr_ind] = ind
        else:
            ind = ind + 1
    return a_seq_pog_dict

def is_subset_extant_dist(an_dict,seq_pog_dict):
    return all(val in seq_pog_dict.get(key, None) for key, val in an_dict.items())


def check_ancestor_indel_distribution(treefile,indelfastafile,extant_alignment_file,add_start_end):
    # load the fasta file
    indel_sequences = FastaFile(indelfastafile)
    # load the tree
    tree_file = open(treefile,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)

    # create extant distribution
    extant_dist_dict = create_extant_distribution(extant_alignment_file)

    # check ancestor indel pattern distribution
    ancestor_in_dist = 0
    all_ancestors = 0

    for n in tree.traverse():
        if n.is_leaf() == False:
            all_ancestors += 1
            new_sequence = indel_sequences.fetch(n.name)
            if add_start_end:
                new_sequence = 'x' + new_sequence + 'x'
            sequence_length = len(new_sequence)
            # check if the sequence is present in the union of all edges
            an_indel_pattern = convert_indel_edge(new_sequence)
            ancestor_in_dist += is_subset_extant_dist(an_indel_pattern,extant_dist_dict)

    #print("Total ancestors in extant distribution", ancestor_in_dist/all_ancestors * 100)
    return round((1 - (ancestor_in_dist/all_ancestors)) * 100,1)


# main function
def main(nwk_file_path,indel_fasta_solution_file,extant_alignment_file,add_start_end):
    # get the parimony score
    out_pattern_percent = check_ancestor_indel_distribution(nwk_file_path,indel_fasta_solution_file,extant_alignment_file,add_start_end)
    print(f"Total Ancestor indel patterns out of extant indel pattern distribution {out_pattern_percent}")

if __name__ == "__main__":
  main(nwk_file_path,indel_fasta_solution_file,extant_alignment_file,add_start_end)

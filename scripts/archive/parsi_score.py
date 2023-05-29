import string
import Levenshtein as lev
from ete3 import Tree
from pysam import FastaFile
import sys

def distance_lev(str1,str2):
    distance_strings = lev.distance(str1,str2)
    return distance_strings

def parsi_score_cal(t):
    # calculate the average distance for each node
    parsi_score = 0
    for n in t.traverse():
        if n.is_leaf() == False:
            neigh_dist_list = []
            # for each neighbour calculate the distance
            for s in n.neighbour_sequence:
                neigh_dist_list.append(distance_lev(n.sequence,s))
            # average distance
            node_avg_distance = sum(neigh_dist_list) / len(neigh_dist_list)
            parsi_score = parsi_score + node_avg_distance
    return parsi_score

def parse_tree_cal_score(nwk_file,ancestor_file,extant_file):
    print("parsing the tree")
    # load tree
    tree_string_file = open(nwk_file,"r")
    my_tree = tree_string_file.read() + ";"
    t = Tree(my_tree, format=1)

    # load the ancestor sequences
    ancestors_sequences_info = FastaFile(ancestor_file)
    # load extant sequences
    extant_sequences_info = FastaFile(extant_file)

    #load sequences on the tree
    for n in t.traverse():
        try:
            n.add_features(sequence = ancestors_sequences_info.fetch(n.name))
        except:
            n.add_features(sequence = extant_sequences_info.fetch(n.name))

    # traverse tree from top down and append 3 sequences ( except for root) for each node to calculate the Lev distance
    # from its neighbours
    print("finding the neighbours")
    for n in t.traverse():
        if n.is_root() == True: # has no parent
            neighbour_sequence_list = []
            for c in n.children:
                neighbour_sequence_list.append(c.sequence)
            n.add_features(neighbour_sequence = neighbour_sequence_list)

        elif n.is_leaf() == False:
            neighbour_sequence_list = []
            for c in n.children:
                neighbour_sequence_list.append(c.sequence)
            neighbour_sequence_list.append(n.up.sequence) # parent sequence
            n.add_features(neighbour_sequence = neighbour_sequence_list)
    print("calculating the score")
    score = parsi_score_cal(t)
    return score


def main():
    # input file and output file
    nwk_file        = sys.argv[1]
    ancestor_file   = sys.argv[2]
    extant_file     = sys.argv[3]
    print("Starting Score Calculation")
    parsimony_score = parse_tree_cal_score(nwk_file,ancestor_file,extant_file)
    print("Tree Parsimony Score ::",parsimony_score)


if __name__ == "__main__":
    main()

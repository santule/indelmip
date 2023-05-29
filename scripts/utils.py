''' Data Processing for MIP model '''

# libraries
from collections import defaultdict
from ete3 import Tree
import numpy as np
from pysam import FastaFile,FastxFile
import pickle
import sys
import pickle

# class with functions to build data structures
class IndelsInfo:
    def __init__(self,fasta_file,nwk_file_path,folder_location):

        self.input_file = fasta_file
        self.nwk_file_path = nwk_file_path
        self.ancestor_list = []
        self.tree_neighbor_dict = defaultdict(list)
        self.ancestor_info = []
        self.sequence_length = 0
        self.folder_location = folder_location
        self.extant_dict = {}
        self.seq_pog_dict = defaultdict(set)
        self.seq_pog_reverse_dict = defaultdict(set)

    # create node types for each position for each sequences
    def create_node_type(self,seq_name):
        node_type_dict = defaultdict(list)
        node_type_dict[(seq_name,'start')] = [0]
        node_type_dict[(seq_name,'end')] = [self.sequence_length - 1]
        all_keys = list(self.seq_pog_dict.keys())
        all_keys.remove(0)
        all_keys.sort()
        node_type_dict[(seq_name,'fwd_back_pos')] = all_keys # remove the start position

        # check if there are any dead nodes
        dead_nodes = list(set(range(1,self.sequence_length - 1)).difference(set(all_keys)))
        node_type_dict[(seq_name,'dead_pos')] = dead_nodes

        return node_type_dict

    # function to find next position that is filled
    def next_pos(self,str1,curr_pos,seq_len):
        start_pos = curr_pos + 1

        while(start_pos < len(str1)):
            if str1[start_pos] != '-':
                return start_pos
            else:
                start_pos = start_pos + 1
        return seq_len

    # function to convert a sequence to adj matrix
    def convert_to_edges(self,seq_str):
        seq_len = len(seq_str)
        ind = 0
        while(ind < seq_len - 1):
            if seq_str[ind] != '-':
                curr_ind = ind
                ind = self.next_pos(seq_str,curr_ind,seq_len - 1) # find the next filled position
                self.seq_pog_dict[curr_ind].add(ind)
                self.seq_pog_reverse_dict[ind].add(curr_ind)
            else:
                ind = ind + 1

    # 1 - convert fasta file to pogs, seq name + binary
    def get_extant_data(self):
        seq_name_list     = []
        seq_binary_list   = []

        with FastxFile(self.input_file) as fh:
            for entry in fh:
                # add start and end string to the sequence
                seq_name = entry.name
                new_sequence = 'x' + entry.sequence + 'x'
                self.sequence_length = len(new_sequence)

                # convert to edges
                self.convert_to_edges(new_sequence)

                # binarise sequences
                seq_binary   = [0 if s == '-' else 1 for s in new_sequence]
                seq_binary   = ''.join([str(x) for x in seq_binary])

                # add to the list
                seq_name_list.append(seq_name)
                seq_binary_list.append(seq_binary)

        self.extant_dict = dict(zip(seq_name_list, seq_binary_list))
        return self.extant_dict

    # 2 - create neighbour dict using the tree file
    def get_tree_data(self):

        ''' create neighbor dict '''
        tree_file = open(self.nwk_file_path,"r")
        my_tree = tree_file.read() + ";"
        tree = Tree(my_tree, format=1)

        # add node names to the internal branches
        # edge = 0
        # for n in tree.traverse():
        #     if not n.is_leaf():
        #         n.name = "NODE_%d" %edge
        #         edge += 1
        #         self.ancestor_list.append(n.name)

        # create neighbourhood object
        for n in tree.traverse():
            if n.is_leaf() == False:
                self.ancestor_list.append(n.name)
                for c in n.children:
                    self.tree_neighbor_dict[n.name] += [c.name]

        return self.tree_neighbor_dict

    # 3 - ancestor data - all ancestors, aggregated pog, aggregated adj mat
    def get_ancestor_data(self):

        # all ancestors name
        ancestor_branchpoints = self.ancestor_list

        # create node type dict
        ancestor_node_type = self.create_node_type('ANCESTOR')
        self.ancestor_info = [ancestor_branchpoints,self.seq_pog_dict,self.seq_pog_reverse_dict,ancestor_node_type]
        return self.ancestor_info

    # save Dataset
    def save_data(self):
        # neighbor dict
        with open(self.folder_location + '/neighbor_dict.pkl','wb') as f:
            pickle.dump(self.tree_neighbor_dict,f)
        # ancestor Info
        with open(self.folder_location + '/ancestor_info.pkl','wb') as f:
            pickle.dump(self.ancestor_info,f)
        # extant info
        with open(self.folder_location + '/extant_data.pkl','wb') as f:
            pickle.dump(self.extant_dict,f)

        return self.extant_dict,self.ancestor_info,self.tree_neighbor_dict

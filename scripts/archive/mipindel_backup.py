''' Create input files for MIP model - input is phylogenetic tree and alignment fasta file '''

# libraries
from collections import defaultdict
from ete3 import Tree
import numpy as np
from pysam import FastaFile,FastxFile
import pickle
import sys
import getopt
import gurobipy as gp
from gurobipy import abs_,quicksum
from gurobipy import GRB
import time
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

# class for mip model building and training
class PhyloTreeMIP:
    def __init__(self,extant_data,ancestor_data,neighbor_dict):

        # Define the configuration  and decision variables for the tree
        self.extant_data = extant_data
        self.ancestor_data = ancestor_data
        self.sequence_length = len(list(self.extant_data.values())[0])
        self.neighbor_dict = neighbor_dict
        self.objective = []
        self.mip_ancestor_fasta_file = 'mip_ancestor_indel.fasta'
        self.all_node_paths = {}
        self.mip_model_file = 'mip_model.mps.zip'
        self.sequence_range = range(1,self.sequence_length-1)
        self.extant_list = self.extant_data.keys()
        self.ancestor_list = self.ancestor_data[0]

        # MIP data structures
        self.ancestorsequence = {}
        self.edges = {}
        self.penalty = {}
        self.diff = {}
        self.objective = []

        # 2 - create a new model
        self.m = gp.Model("PreferredPathSolve")

    def add_pos_constraints_ancestors(self):
        ''' function to add constraints for positions in ancestors  '''

        # V - create variable for each ancestor sequence
        for ancestor in self.ancestor_list:
            an = self.m.addVars(self.sequence_range,vtype=GRB.BINARY)
            self.ancestorsequence[ancestor] = an

    def add_edge_constraints_ancestors(self):
        ''' function to add constraints for edges in ancestors  '''

        ancestor_fwd_edges  = self.ancestor_data[1]
        ancestor_bkwd_edges = self.ancestor_data[2]
        ancestor_node_type  = self.ancestor_data[3]

        # for each ancestor
        for ancestor in self.ancestor_list:

            # START NODES
            for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','start')):

                all_edges_from_pos = []
                for pos_to in ancestor_fwd_edges[fwd_back_pos]:
                    edge_id = (ancestor,fwd_back_pos,pos_to)
                    # V - var for each edge from start node
                    e = self.m.addVar(vtype=GRB.BINARY)
                    self.edges[edge_id] = e
                    all_edges_from_pos.append(e)

                # C - 1 edge has to be can be used
                self.m.addConstr(quicksum(all_edges_from_pos) == 1)

            # DEAD NODES
            if ancestor_node_type.get(('ANCESTOR','dead_pos')) :
                for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','dead_pos')):
                    pos_id = (ancestor,fwd_back_pos)
                    # C - constraint for dead pos to 0 as they cannot find complete path
                    self.m.addConstr(self.ancestorsequence[ancestor][fwd_back_pos] == 0)

            # FULLY CONNECTED NODES
            for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','fwd_back_pos')):

                pos_id = (ancestor,fwd_back_pos)
                # C - edges going out of the node
                all_edges_from_pos = []
                for pos_to in ancestor_fwd_edges[fwd_back_pos]:
                    edge_id = (ancestor,fwd_back_pos,pos_to)
                    # V - var for each edge going from the node
                    e = self.m.addVar(vtype=GRB.BINARY)
                    self.edges[edge_id] = e
                    all_edges_from_pos.append(e)

                # C - edges coming in into the node
                edges_coming_in_list = []
                for edges_coming_in_item in ancestor_bkwd_edges[fwd_back_pos]:
                    edge_to_id = (ancestor,edges_coming_in_item,fwd_back_pos)
                    edges_coming_in_list.append(self.edges[edge_to_id])


                # C - sum(edges) = position
                self.m.addConstr(quicksum(all_edges_from_pos)   == self.ancestorsequence[ancestor][fwd_back_pos])
                # C - sum(edges) = position
                self.m.addConstr(quicksum(edges_coming_in_list) == self.ancestorsequence[ancestor][fwd_back_pos])

            # END NODES
            for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','end')):

                all_edges_to_pos = []
                for pos_from in ancestor_bkwd_edges[fwd_back_pos]:
                    edge_id = (ancestor,pos_from,fwd_back_pos)
                    all_edges_to_pos.append(self.edges[edge_id])

                # C - sum(edges) = position
                self.m.addConstr(quicksum(all_edges_to_pos) == 1)

    def add_penalty_constraint(self):
        ''' function to add penalty constraints for whole tree  '''

        # difference constraints
        for node,node_neighbor in self.neighbor_dict.items():
            for node_neighbor_item in node_neighbor:

                # V - add position difference variable
                node_pos_var          = self.ancestorsequence[node]
                nn_pair               = (node,node_neighbor_item)
                pen                   = self.m.addVars(self.sequence_range,vtype=GRB.BINARY)
                self.penalty[nn_pair] = pen

                # check if it is a extant node or not
                if node_neighbor_item in self.extant_list:
                    node_neighbor_pos_var = [int(i) for i in self.extant_data[node_neighbor_item]]
                else:
                    node_neighbor_pos_var = self.ancestorsequence[node_neighbor_item]
                    diff_pos              = self.m.addVars(self.sequence_range,vtype=GRB.BINARY)
                    self.diff[nn_pair]    = diff_pos


                for pos in range(1,self.sequence_length - 1):

                    if node_neighbor_item in self.extant_list:
                        # C - difference varaibles only if it is not an extant
                        if node_neighbor_pos_var[pos] == 1:
                            self.diff[(node,node_neighbor_item,pos)]  = 1 - node_pos_var[pos]
                        elif node_neighbor_pos_var[pos] == 0:
                            self.diff[(node,node_neighbor_item,pos)]  = node_pos_var[pos]


                        # C - penalty constraits
                        if pos == 1: # penalty for first position is simple
                            self.m.addConstr(pen[pos] == self.diff[(node,node_neighbor_item,pos)])
                        else:
                            self.m.addConstr(pen[pos] >= self.diff[(node,node_neighbor_item,pos)] -\
                                             self.diff[(node,node_neighbor_item,pos-1)])

                        # O - add difference to the objective
                        self.objective.append(self.diff[(node,node_neighbor_item,pos)])

                        # O - add penalty to the objective
                        self.objective.append(2 * pen[pos])

                    else:
                        # C - Abs constraints
                        self.m.addConstr( diff_pos[pos] <= node_pos_var[pos] + node_neighbor_pos_var[pos])
                        self.m.addConstr( diff_pos[pos] >= node_pos_var[pos] - node_neighbor_pos_var[pos])
                        self.m.addConstr( diff_pos[pos] >= node_neighbor_pos_var[pos] - node_pos_var[pos])
                        self.m.addConstr( diff_pos[pos] <= 2 - node_neighbor_pos_var[pos] - node_pos_var[pos])

                        # C - penalty constraits
                        if pos == 1: # penalty for first position is simple
                            self.m.addConstr(pen[pos] == diff_pos[pos])
                        else:
                            self.m.addConstr(pen[pos] >= diff_pos[pos] - diff_pos[pos-1])

                        # O - add difference to the objective
                        self.objective.append(diff_pos[pos])

                        # O - add penalty to the objective
                        self.objective.append(2 * pen[pos])

    def train(self,n_threads,time_out):
        # Params
        self.m.Params.Threads = n_threads
        self.m.Params.TimeLimit = time_out*60
        self.m.Params.LogToConsole = 0
        self.m.Params.Method = 1
        #self.m.Params.LogFile = "mip_model.log"

        # Optimize
        self.total_objective = quicksum(self.objective)
        self.m.setObjective(self.total_objective, GRB.MINIMIZE)
        self.m.update()

        # save the model
        #self.m.write(self.mip_model_file )
        self.m.optimize()

        #Is feasible?
        return self.m.SolCount > 0

    def get_info(self):
        info_all = {}
        info_all["objective"] = self.m.ObjVal
        info_all["bound"] = self.m.ObjBound
        info_all["gap"] = self.m.MIPGap
        info_all["is_optimal"] = (self.m.status == GRB.OPTIMAL)
        info_all["num_nodes"] = self.m.NodeCount
        info_all["num_vars"] = self.m.numVars

        if self.m.SolCount > 0:
            print("objective: %0.2f"%info_all["objective"])
            print("bound: %0.2f"%info_all["bound"])
            print("gap: %0.2f"%info_all["gap"])

        return info_all

    def output_fasta(self,all_node_paths):
        # convert output file to FASTA file
        with open(self.mip_ancestor_fasta_file,mode='w') as fout:
            for node_name,sequence in all_node_paths.items():
                fout.write('>' + str(node_name) + '\n')
                sequence_str = ''.join([str(s) for s in sequence])
                fout.write(str(sequence_str) + '\n')

    def get_solution(self):
        # get the path for extants - should be same as the input
        all_node_paths = {}
        for sequence_name,preferred_path in self.extant_data.items():
            all_node_paths[sequence_name] = preferred_path

        # get the path for ancestor
        for ancestor in self.ancestor_data[0]:
            preferred_path = []
            preferred_path.append('1') # start position
            for pos in range(1,self.sequence_length - 1): # remove the start and end fake positions
                preferred_path.append(str(int(self.ancestorsequence[ancestor][pos].X)))
            preferred_path.append('1') # end position
            all_node_paths[ancestor] = "".join(preferred_path)
        return all_node_paths

# process input data for MIP
def data_processing(extant_sequence_file,nwk_file_path,folder_location):

    # prepare input Dataset
    print("Processing Input Files")
    MIPIndel      = IndelsInfo(extant_sequence_file,nwk_file_path,folder_location) # class
    print("1 - Processing Extant Data")
    extant_data   = MIPIndel.get_extant_data() # pog, adj matrix, binary, node type, sequence name
    print("2 - Preparing Tree Data")
    neighbor_dict = MIPIndel.get_tree_data() # neighbor info
    print("3 - Preparing Ancestor Data")
    ancestor_data = MIPIndel.get_ancestor_data() # ancestor list, ancestor pog, node type
    print("4 - Saving Data")
    extant_data,ancestor_data,neighbor_dict = MIPIndel.save_data() # save data
    print("Done")

    # Info about the data
    total_sequences = len(ancestor_data[0]) + 1
    print("Total extant sequences",total_sequences)
    print("Sequence length",len(list(extant_data.values())[0]))
    return extant_data,ancestor_data,neighbor_dict

# process mip
def mip_processing(extant_data,ancestor_data,neighbor_dict):
    # Create the MIP model
    start = time.time()
    print("Start Time:",start)
    # initialise the class
    PyTree = PhyloTreeMIP(extant_data,ancestor_data,neighbor_dict)
    # variable position constraints for ancestors
    print("Adding Ancestor Constraints")
    PyTree.add_pos_constraints_ancestors()
    # edge constraints for ancestors
    print("Adding Edges Constraints")
    PyTree.add_edge_constraints_ancestors()
    # position difference constraints
    print("Adding Penalty Constraints")
    PyTree.add_penalty_constraint()
    total_time = ((time.time()-start)/60)
    print("-----------------------------")
    print("Total time to create model = %0.2f[mins]"%total_time)

    # Run MIP
    print('Training MIP Model')
    start = time.time()
    print("Start Time:",start)
    n_threads = 1
    time_out = 60 # 1 hour

    is_sat = PyTree.train(n_threads, time_out)
    print("is_sat",is_sat)
    total_time = ((time.time()-start)/60)
    print("-----------------------------")
    print("Total time to solve model= %0.2f[mins]"%total_time)
    info = PyTree.get_info()
    info["total_time"] = total_time
    info["is_sat"] = is_sat
    print("info",info)

    if is_sat:
        all_node_paths = PyTree.get_solution()
        PyTree.output_fasta(all_node_paths)
        print(f"MIP indel output file created and saved")
    else:
        print("Did not find any satisfactory solution to the model")

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('mipindel.py -a alignment file | -n newick tree file | -o output folder')
    exit()

# main function
def main():

    extant_sequence_file = None
    nwk_file_path = None
    output_folder_location = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:n:o:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-a']:
            extant_sequence_file = arg
        elif opt in ['-n']:
            nwk_file_path = arg
        elif opt in ['-o']:
            output_folder_location = arg

    if extant_sequence_file is None or nwk_file_path is None or output_folder_location is None:
        help()

    # 1 - process the input data for MIP
    extant_data,ancestor_data,neighbor_dict = data_processing(extant_sequence_file,nwk_file_path,output_folder_location)
    # 2 - build model
    mip_processing(extant_data,ancestor_data,neighbor_dict)


if __name__ == "__main__":
  main()

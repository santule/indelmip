''' Run MIP for indel inference '''

# libraries
import gurobipy as gp
from gurobipy import abs_,quicksum
from gurobipy import GRB
import time
import json
from collections import defaultdict
from ete3 import Tree
import numpy as np
from pysam import FastaFile,FastxFile
import re
import torch
from torch.utils.data import Dataset



# pytorch dataset to save the extant data
class AjMat_Dataset(Dataset):
    def __init__(self,adj_mat,seq_name,seq_fwd_pog,seq_rvs_pog,node_type,seq_binary):
        self.adj_mat = adj_mat
        self.seq_name = seq_name
        self.seq_fwd_pog = seq_fwd_pog
        self.seq_rvs_pog = seq_rvs_pog
        self.node_type = node_type
        self.seq_binary = seq_binary

    def __len__(self):
        return len(self.adj_mat)
    def __getitem__(self, idx):
        return self.adj_mat[idx],self.seq_name[idx],self.seq_fwd_pog[idx],\
               self.seq_rvs_pog[idx],self.node_type[idx],self.seq_binary[idx]

class IndelsInfo:
    def __init__(self,AjMat_Dataset,fasta_file,nwk_file_path):

        self.AjMat_Dataset = AjMat_Dataset
        self.input_file = fasta_file
        self.nwk_file_path = nwk_file_path
        self.ancestor_list = []
        self.tree_neighbor_dict = defaultdict(list)
        self.ancestor_info = []
        self.sequence_length = 0
        self.Extant_AdjMat_dataset = AjMat_Dataset


    # create node types for each position for each sequences
    def create_node_type(self, seq_fwd_pog,seq_rvs_pog,seq_name):
        node_type_dict = defaultdict(list)
        node_type_dict[(seq_name,'start')] = [0]
        node_type_dict[(seq_name,'end')] = [self.sequence_length - 1]

        for n in range(1,self.sequence_length - 1):
            if n in seq_fwd_pog.keys() and n in seq_rvs_pog.keys(): #if node has forward and backward
                node_type_dict[(seq_name,'fwd_back_pos')] += [n]
            elif n in seq_fwd_pog.keys() and n not in seq_rvs_pog.keys():
                node_type_dict[(seq_name,'fwd_pos')] += [n]
            elif n not in seq_fwd_pog.keys() and n in seq_rvs_pog.keys():
                node_type_dict[(seq_name,'back_pos')] += [n]
            else:
                node_type_dict[(seq_name,'dead_pos')] += [n]
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
    def convert_to_adj_mat(self,seq_str):
        seq_len = len(seq_str)
        aj_mat_array = np.zeros((seq_len,seq_len))

        next_filled = []
        ind = 0

        while(ind < seq_len - 1):
            if seq_str[ind] != '-':
                curr_ind = ind
                ind = self.next_pos(seq_str,curr_ind,seq_len - 1) # find the next filled position
                next_filled.append((curr_ind,ind))
                aj_mat_array[curr_ind,ind] = 1
            else:
                ind = ind + 1
        return aj_mat_array

    # convert adj matrix into pog dictionary
    def create_extant_pog(self,adj_mat_t):
        x_summ = np.column_stack(np.where(adj_mat_t))
        seq_fwd_pog_dict = dict(zip(x_summ[:,0], x_summ[:,1]))
        seq_rvs_pog_dict = dict(zip(x_summ[:,1], x_summ[:,0]))
        return seq_fwd_pog_dict,seq_rvs_pog_dict

    # 1 - convert fasta file to adj matrix, pog, node type, seq binary into pytorch dataset
    def get_extant_data(self):
        adj_mat_list      = []
        seq_name_list     = []
        seq_fwd_pog_list  = []
        node_type_list    = []
        seq_rev_pog_list  = []
        seq_binary_list   = []

        with FastxFile(self.input_file) as fh:
            for entry in fh:
                # add start and end string to the sequence
                seq_name = entry.name
                new_sequence = 'x' + entry.sequence + 'x'
                self.sequence_length = len(new_sequence)

                # convert to adj matrix
                seq_adj_mat  = self.convert_to_adj_mat(new_sequence)

                # binarise sequences
                seq_binary   = ''.join(sum(seq_adj_mat).astype(int).astype(str))
                # make start pos as 1 for start node
                seq_binary = '1' + seq_binary[1:]

                # convert to pog structure
                seq_fwd_pog,seq_rvs_pog = self.create_extant_pog(seq_adj_mat)

                # create node type dict
                node_type = self.create_node_type(seq_fwd_pog,seq_rvs_pog,seq_name)

                # add to the list
                #adj_mat_t = torch.from_numpy(seq_adj_mat)
                adj_mat_list.append(seq_adj_mat)
                seq_name_list.append(seq_name)
                seq_fwd_pog_list.append(seq_fwd_pog)
                seq_rev_pog_list.append(seq_rvs_pog)
                node_type_list.append(node_type)
                seq_binary_list.append(seq_binary)

        # save it into pytorch dataset
        self.Extant_AdjMat_dataset = self.AjMat_Dataset(adj_mat_list, seq_name_list, seq_fwd_pog_list,
                                              seq_rev_pog_list, node_type_list,seq_binary_list)

        return self.Extant_AdjMat_dataset

    # 2 - create neighbour dict using the tree file
    def get_tree_data(self):

        ''' create neighbor dict '''
        tree_file = open(self.nwk_file_path,"r")
        my_tree = tree_file.read() + ";"
        tree = Tree(my_tree, format=1)

        # add node names to the internal branches
        edge = 0
        for n in tree.traverse():
            if not n.is_leaf():
                n.name = "NODE_%d" %edge
                edge += 1
                self.ancestor_list.append(n.name)

        # create neighbourhood object
        for n in tree.traverse():
            if n.is_leaf() == False:
                for c in n.children:
                    self.tree_neighbor_dict[n.name] += [c.name]

        return self.tree_neighbor_dict

    # 3 - ancestor data - all ancestors, aggregated pog, aggregated adj mat
    def get_ancestor_data(self):

        # all ancestors name
        ancestor_branchpoints = self.ancestor_list

        ancestor_fwd_pog = defaultdict(list)
        ancestor_rvs_pog = defaultdict(list)

        # ancestor adj mat
        ancestor_adj_mat = np.where(sum(self.Extant_AdjMat_dataset[:][0]))

        # ancestor foward and backward pog
        row_col_sum = np.column_stack(np.where(sum(self.Extant_AdjMat_dataset[:][0])))
        for r in row_col_sum:
            pos = r[0]
            next_pos = r[1]
            ancestor_fwd_pog[pos] += [next_pos]
            ancestor_rvs_pog[next_pos] += [pos]

        # create node type dict
        ancestor_node_type = self.create_node_type(ancestor_fwd_pog,ancestor_rvs_pog,'ANCESTOR')
        self.ancestor_info = [ancestor_branchpoints,ancestor_fwd_pog,ancestor_rvs_pog,ancestor_node_type]
        return self.ancestor_info

class PhyloTreeMIP:
    def __init__(self,extant_data,ancestor_data,tree_name,neighbor_dict,mip_ancestor_fasta_file):

        # Define the configuration  and decision variables for the tree
        self.extant_data = extant_data
        self.ancestor_data = ancestor_data
        self.tree_name = tree_name
        self.sequence_length = len(self.extant_data[0][5])
        self.neighbor_dict = neighbor_dict
        self.objective = []
        self.mip_ancestor_fasta_file = mip_ancestor_fasta_file
        self.all_node_paths = {}

        # MIP data structures
        self.edges = {}
        self.positions = {}
        self.penalty = {}
        self.diff = {}
        self.objective = []
        self.M = 999

        # 2 - create a new model
        self.m = gp.Model("PreferredPathSolve")


    def add_pos_constraints_extants(self):
        ''' function to add constraints for positions in extants to fix them '''

        for extant_info in self.extant_data:
            sequence_binary = extant_info[5]
            sequence_name   = extant_info[1]
            # V - create variable for each position for each extant sequence
            for pos_idx in range(0,len(sequence_binary)):
                pos_id = (sequence_name,pos_idx)
                pos = self.m.addVar(vtype=GRB.BINARY, name="p-%s-%s"%pos_id)
                self.positions[pos_id] = pos

                # C - fix the extant positions as per binary sequence
                self.m.addConstr(self.positions[pos_id] == int(sequence_binary[pos_idx]),\
                                     name="extant_position_constraint-%s-%s"%pos_id)

    def add_pos_constraints_ancestors(self):
        ''' function to add constraints for positions in ancestors  '''

        ancestor_list = self.ancestor_data[0]
        # V - create variable for each position for each ancestor sequence
        for ancestor in ancestor_list:
            for pos_idx in range(0,self.sequence_length):
                pos_id = (ancestor,pos_idx)
                pos = self.m.addVar(vtype=GRB.BINARY, name="p-%s-%s"%pos_id)
                self.positions[pos_id] = pos

            # C - start pos is always 1
            pos_id = (ancestor,0)
            pos = self.m.addVar(vtype=GRB.BINARY, name="p-%s-%s"%pos_id)
            self.positions[pos_id] = pos
            self.m.addConstr(self.positions[pos_id] == 1,\
                                     name="ancestor_start_end_position_constraint-%s-%s"%pos_id)

            # C - end pos is always 1
            pos_id = (ancestor,self.sequence_length-1)
            pos = self.m.addVar(vtype=GRB.BINARY, name="p-%s-%s"%pos_id)
            self.positions[pos_id] = pos
            self.m.addConstr(self.positions[pos_id] == 1,\
                                     name="ancestor_start_end_position_constraint-%s-%s"%pos_id)

    def add_edge_constraints_ancestors(self):
        ''' function to add constraints for edges in ancestors  '''

        ancestor_list = self.ancestor_data[0]
        ancestor_node_type = self.ancestor_data[3]
        ancestor_fwd_edges = self.ancestor_data[1]
        ancestor_bkwd_edges = self.ancestor_data[2]

        # constraint for each ancestor
        for ancestor in ancestor_list:

            # START NODES
            for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','start')):

                all_edges_from_pos = []
                for pos_to in ancestor_fwd_edges[fwd_back_pos]:
                    edge_id = (ancestor,fwd_back_pos,pos_to)
                    # V - var for each edge from start node
                    e = self.m.addVar(vtype=GRB.BINARY, name='e-%s-%s-%s'%edge_id)
                    self.edges[edge_id] = e
                    all_edges_from_pos.append(e)

                # C - only 1 edge can be used
                pos_id = (ancestor,fwd_back_pos)
                self.m.addConstr(quicksum(all_edges_from_pos) <= 1,name=\
                                          "ancestor_start_edge_constraint-%s-%s"%pos_id)
                # C - sum(edges) = position
                self.m.addConstr(quicksum(all_edges_from_pos) == self.positions[pos_id],\
                                            name = "ancestor_edge_node_recon_constraint-%s-%s"%pos_id)


            # END NODES
            for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','end')):

                all_edges_to_pos = []
                for pos_from in ancestor_bkwd_edges[fwd_back_pos]:
                    edge_id = (ancestor,pos_from,fwd_back_pos)
                    # V - var for each edge to end node
                    e = self.m.addVar(vtype=GRB.BINARY, name='e-%s-%s-%s'%edge_id)
                    self.edges[edge_id] = e
                    all_edges_to_pos.append(e)

                # C - only 1 edge can be used
                pos_id = (ancestor,fwd_back_pos)
                self.m.addConstr(quicksum(all_edges_to_pos) <= 1,name=\
                                          "ancestor_end_edge_constraint-%s-%s"%pos_id)
                # C - sum(edges) = position
                self.m.addConstr(quicksum(all_edges_to_pos) == self.positions[pos_id],\
                                            name = "ancestor_edge_node_recon_constraint-%s-%s"%pos_id)


            # DEAD NODES
            if ancestor_node_type.get(('ANCESTOR','dead_pos')) :
                for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','dead_pos')):
                    pos_id = (ancestor,fwd_back_pos)
                    # C - constraint for dead pos to 0 as they cannot find complete path
                    self.m.addConstr(self.positions[pos_id] == 0,\
                                                name = "ancestor_dead_pos-%s-%s"%pos_id)

            # FULLY CONNECTED NODES
            for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','fwd_back_pos')):


                all_edges_from_pos = []
                for pos_to in ancestor_fwd_edges[fwd_back_pos]:
                    edge_id = (ancestor,fwd_back_pos,pos_to)
                    # V - var for each edge going from the node
                    e = self.m.addVar(vtype=GRB.BINARY, name='e-%s-%s-%s'%edge_id)
                    self.edges[edge_id] = e
                    all_edges_from_pos.append(e)

                # C - only 1 edge can be used
                pos_id = (ancestor,fwd_back_pos)
                self.m.addConstr(quicksum(all_edges_from_pos) <= 1,name=\
                                          "ancestor_edge_constraint-%s-%s"%pos_id)

                # C - sum(edges going in) = sum(edges going out) for each node
                edges_coming_in_list = []
                for edges_coming_in_item in ancestor_bkwd_edges[fwd_back_pos]:
                    edge_to_id = (ancestor,edges_coming_in_item,fwd_back_pos)
                    edges_coming_in_list.append(self.edges[edge_to_id])

                self.m.addConstr(quicksum(edges_coming_in_list) == quicksum(all_edges_from_pos),\
                                                 name="ancestor_edge_recon_constraint-%s-%s"%pos_id)
                # C - sum(edges) = position
                self.m.addConstr(quicksum(all_edges_from_pos) == self.positions[pos_id],\
                                            name="ancestor_edge_node_recon_constraint1-%s-%s"%pos_id)
                # C - sum(edges) = position
                self.m.addConstr(quicksum(edges_coming_in_list) == self.positions[pos_id],\
                                            name="ancestor_edge_node_recon_constraint2-%s-%s"%pos_id)


    # penalty constraint
    def add_penalty_constraint(self):
        ''' function to add penalty constraints for whole tree  '''

        # difference constraints
        for node,node_neighbor in self.neighbor_dict.items():
            for node_neighbor_item in node_neighbor:
                for pos in range(1,self.sequence_length - 1): # penalty start from 1st position only

                    # V - penalty variables for node to node for each position
                    pen_id = (node,node_neighbor_item,pos)
                    pen = self.m.addVar(vtype=GRB.BINARY, name='pe-%s-%s-%s'%pen_id)
                    self.penalty[pen_id] = pen

                    # V - add position difference variable
                    node_pos_var = self.positions[(node,pos)]
                    node_neighbor_pos_var = self.positions[(node_neighbor_item,pos)]
                    diff_id = (node,node_neighbor_item,pos)
                    diff_pos = self.m.addVar(vtype=GRB.BINARY, name='d-%s-%s-%s'%diff_id)
                    self.diff[diff_id] = diff_pos

                    # C - abs difference constraint
                    self.m.addConstr( diff_pos <= node_pos_var + node_neighbor_pos_var,name=\
                                     "diff_constraint_1-%s-%s-%s"%(node,node_neighbor_item,pos))
                    self.m.addConstr( diff_pos >= node_pos_var - node_neighbor_pos_var,name=\
                                     "diff_constraint_2-%s-%s-%s"%(node,node_neighbor_item,pos))
                    self.m.addConstr( diff_pos >= node_neighbor_pos_var - node_pos_var,name=\
                                     "diff_constraint_3-%s-%s-%s"%(node,node_neighbor_item,pos))
                    self.m.addConstr( diff_pos <= 2 - node_neighbor_pos_var - node_pos_var,name=\
                                     "diff_constraint_4-%s-%s-%s"%(node,node_neighbor_item,pos))
                    # O - add objective
                    self.objective.append(diff_pos)

        # gap penalty constraints
        for node,node_neighbor in self.neighbor_dict.items():
            for node_neighbor_item in node_neighbor:
                for pos in range(1,self.sequence_length - 1):  # no penalty for start and end
                    diff_id  = (node,node_neighbor_item,pos)
                    pen_id   = (node,node_neighbor_item,pos)
                    diff_var = self.diff[diff_id]
                    pen_var  = self.penalty[pen_id]

                    if pos == 1: # penalty for first position is simple
                        self.m.addConstr(pen_var == diff_var,"penalty_constraint-%s-%s-%s"%\
                                         (node,node_neighbor_item,pos))
                    else:
                        pen_prev_id = (node,node_neighbor_item,pos - 1)
                        prev_pen_var =  self.penalty[pen_prev_id]
                        prev_diff_var = self.diff[pen_prev_id]

                        # C - gap opening penalty
                        self.m.addConstr(diff_var - prev_diff_var >= 1 - self.M * (1 - pen_var),\
                                         name="penalty_constraint_1-%s-%s-%s"%\
                                         (node,node_neighbor_item,pos))
                        self.m.addConstr(diff_var - prev_diff_var <= self.M * (pen_var),\
                                         name="penalty_constraint_2-%s-%s-%s"%\
                                         (node,node_neighbor_item,pos))

                    # O - add penalty to the objective
                    self.objective.append(2 * pen_var)



    def train(self,n_threads,time_out):
        # Params
        self.m.Params.Threads = n_threads
        self.m.Params.TimeLimit = time_out*60
        self.m.Params.LogToConsole = 0
        self.m.Params.Degenmoves=0

        # Optimize
        self.total_objective = sum([o for o in self.objective])
        self.m.setObjective(self.total_objective, GRB.MINIMIZE)
        self.m.update()

        self.m.write(('mip_model_' + self.tree_name + '.lp'))
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
        info_all["num_vars"] = self.m.NumIntVars + self.m.NumBinVars

        if self.m.SolCount > 0:
            print("objective: %0.2f"%info_all["objective"])
            print("bound: %0.2f"%info_all["bound"])
            print("gap: %0.2f"%info_all["gap"])

        return info_all

    def get_solution(self):

        # get the path for extants - should be same as the input
        for extant_info in self.extant_data:
            sequence_name   = extant_info[1]
            preferred_path = []
            for pos in range(0,self.sequence_length):
                pos_id = (sequence_name,pos)
                preferred_path.append(int(self.positions[pos_id].X))
            self.all_node_paths[sequence_name] = preferred_path

        # get the path for ancestor
        for ancestor in self.ancestor_data[0]:
            preferred_path = []
            for pos in range(0,self.sequence_length):
                pos_id = (ancestor,pos)
                preferred_path.append(int(self.positions[pos_id].X))
            self.all_node_paths[ancestor] = preferred_path

        # get the differnece and penalty solution
        score_dict = {}
        overall_score = 0
        for node,node_neighbor in self.neighbor_dict.items():
            for node_neighbor_item in node_neighbor:
                total_score = 0
                for pos in range(1,self.sequence_length - 1): # penalty start from 1st position only
                    pen_id  = (node,node_neighbor_item,pos)
                    diff_id = (node,node_neighbor_item,pos)

                    total_score = total_score + 2 * int(self.penalty[pen_id].X)
                    total_score = total_score + int(self.diff[diff_id].X)

                score_dict[(node,node_neighbor_item)] = total_score
                overall_score = overall_score + total_score
        return self.all_node_paths,score_dict

    def output_fasta(self):
        # convert output file to FASTA file
        with open(self.mip_ancestor_fasta_file,mode='w') as fout:
            for node_name,sequence in self.all_node_paths.items():
                fout.write('>' + str(node_name) + '\n')
                sequence_str = ''.join([str(s) for s in sequence])
                fout.write(str(sequence_str) + '\n')


def main():
    folder_location         = '/Users/sanjanatule/Documents/uq/Projects/MIPIndel/scripts/mip_files/'


    ## Sample tree 1
    nwk_file_path           = '/media/WorkingSpace/Share/mipindel/data/st1/input_tree.nwk'
    extant_sequence_file    = '/media/WorkingSpace/Share/mipindel/data/st1/input_extants.fasta'
    mip_ancestor_fasta_file = '/media/WorkingSpace/Share/mipindel/data/st1/mip_ancestor_indel.fasta'
    tree_name               = 'sample1'

    ## CYP2U - 165
    # nwk_file_path           = '/media/WorkingSpace/Share/mipindel/data/CYP2U_165/CYP2U_165.nwk'
    # extant_sequence_file    = '/media/WorkingSpace/Share/mipindel/data/CYP2U_165/CYP2U_165.aln'
    # mip_ancestor_fasta_file = "/media/WorkingSpace/Share/mipindel/data/CYP2U_165/mip_ancestor_indel.fasta"
    # tree_name = 'cyp2u_165'

    # ## CYP2U - 359
    # nwk_file_path           = '/media/WorkingSpace/Share/mipindel/data/CYP2U_359/CYP2U_359.nwk'
    # extant_sequence_file    = '/media/WorkingSpace/Share/mipindel/data/CYP2U_359/CYP2U_359.aln'
    # mip_ancestor_fasta_file = "/media/WorkingSpace/Share/mipindel/data/CYP2U_359/mip_ancestor_indel.fasta"
    # tree_name = 'cyp2u_359'

    # ## DHAD - 1612
    # nwk_file_path           = '/media/WorkingSpace/Share/mipindel/data/DHAD_1612/DHAD_1612.nwk'
    # extant_sequence_file    = '/media/WorkingSpace/Share/mipindel/data/DHAD_1612/DHAD_1612.aln'
    # mip_ancestor_fasta_file = "/media/WorkingSpace/Share/mipindel/data/DHAD_1612/mip_ancestor_indel.fasta"
    # tree_name = 'DHAD_1612'

    ## MBL
    # nwk_file_path           = '/media/WorkingSpace/Share/mipindel/data/MBL/nuclease_filt_i10.aln.treefile.nwk'
    # extant_sequence_file    = '/media/WorkingSpace/Share/mipindel/data/MBL/nuclease_filt_i10.aln'
    # mip_ancestor_fasta_file = "/media/WorkingSpace/Share/mipindel/data/MBL/mip_ancestor_indel.fasta"
    # tree_name = 'MBL'

    # CYPU - Anthony
    # nwk_file_path           = '/media/WorkingSpace/Share/mipindel/data/anthony/CYP19_Putative_6_DASH.nwk'
    # extant_sequence_file    = '/media/WorkingSpace/Share/mipindel/data/anthony/CYP19_Putative_6_DASH.fasta'
    # mip_ancestor_fasta_file = "/media/WorkingSpace/Share/mipindel/data/anthony/mip_ancestor_indel.fasta"
    # tree_name = 'anthony'

    # prepare input Dataset
    MIPIndel      = IndelsInfo(AjMat_Dataset,extant_sequence_file,nwk_file_path) # class
    extant_data   = MIPIndel.get_extant_data() # pog, adj matrix, binary, node type, sequence name
    neighbor_dict = MIPIndel.get_tree_data() # neighbor info
    ancestor_data = MIPIndel.get_ancestor_data() # ancestor list, ancestor pog, node type

    # Info about the data
    total_sequences = len(ancestor_data[0]) + 1
    print("TOTAL EXTANT SEQUENCES",total_sequences)
    print("SEQUENCE LENGTH",len(extant_data[0][5]))

    # Run MIP
    start = time.time()
    print("Start Time:",start)
    n_threads = 1
    time_out = 60 # 1 hour

    # initialise the class
    PyTree = PhyloTreeMIP(extant_data,ancestor_data,tree_name,neighbor_dict,mip_ancestor_fasta_file)

    # variable and position constraints for extants
    PyTree.add_pos_constraints_extants()
    # variable and position constraints for ancestors
    PyTree.add_pos_constraints_ancestors()
    # variable and edge constraints for ancestors
    PyTree.add_edge_constraints_ancestors()
    # position difference constraints
    PyTree.add_penalty_constraint()

    is_sat = PyTree.train(n_threads, time_out)
    print("is_sat",is_sat)
    total_time = ((time.time()-start))
    print("-----------------------------")
    print("Total time = %0.2f[m]"%total_time)
    info = PyTree.get_info()
    info["total_time"] = total_time
    info["is_sat"] = is_sat
    print("info",info)

    if is_sat:
        all_node_paths,score_dict = PyTree.get_solution()
        PyTree.output_fasta()
    else:
        print("Did not find any satisfactory solution to the model")

if __name__ == "__main__":
  main()

''' Run MIP for indel inference '''

# libraries
import gurobipy as gp
from gurobipy import abs_,quicksum
from gurobipy import GRB
import time
import json
from collections import defaultdict
import numpy as np
import torch
from torch.utils.data import Dataset
import pickle



# pytorch dataset to save the extant data
class AjMat_lean_Dataset(Dataset):
    def __init__(self,seq_name,seq_binary):
        self.seq_name = seq_name
        self.seq_binary = seq_binary

    def __len__(self):
        return len(self.seq_name)
    def __getitem__(self, idx):
        return self.seq_name[idx],self.seq_binary[idx]


class PhyloTreeMIP:
    def __init__(self,extant_data,ancestor_data,tree_name,neighbor_dict,mip_ancestor_fasta_file,mip_model_file):

        # Define the configuration  and decision variables for the tree
        self.extant_data = extant_data
        self.ancestor_data = ancestor_data
        self.tree_name = tree_name
        self.sequence_length = len(list(self.extant_data.values())[0])
        self.neighbor_dict = neighbor_dict
        self.objective = []
        self.mip_ancestor_fasta_file = mip_ancestor_fasta_file
        self.all_node_paths = {}
        self.mip_model_file = mip_model_file

        # MIP data structures
        self.edges = {}
        self.positions = {}
        self.penalty = {}
        self.diff = {}
        self.objective = []
        self.M = 999

        # 2 - create a new model
        self.m = gp.Model("PreferredPathSolve")


    # def add_pos_constraints_extants(self):
    #     ''' function to add constraints for positions in extants to fix them '''
    #
    #     for extant_info in self.extant_data:
    #         sequence_binary = extant_info[1]
    #         sequence_name   = extant_info[0]
    #         # V - create variable for each position for each extant sequence
    #         for pos_idx in range(0,len(sequence_binary)):
    #             pos_id = (sequence_name,pos_idx)
    #             pos = self.m.addVar(vtype=GRB.BINARY, name="p-%s-%s"%pos_id)
    #             self.positions[pos_id] = pos
    #
    #             # C - fix the extant positions as per binary sequence
    #             self.m.addConstr(self.positions[pos_id] == int(sequence_binary[pos_idx]),\
    #                                  name="extant_position_constraint-%s-%s"%pos_id)

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
            self.m.addConstr(self.positions[pos_id] == 1,\
                                     name="ancestor_start_end_position_constraint-%s-%s"%pos_id)

            # C - end pos is always 1
            pos_id = (ancestor,self.sequence_length-1)
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

            # END NODES
            for fwd_back_pos in ancestor_node_type.get(('ANCESTOR','end')):

                 all_edges_to_pos = []
                 for pos_from in ancestor_bkwd_edges[fwd_back_pos]:
                     edge_id = (ancestor,pos_from,fwd_back_pos)
                     all_edges_to_pos.append(self.edges[edge_id])

                 # C - only 1 edge can be used
                 pos_id = (ancestor,fwd_back_pos)
                 self.m.addConstr(quicksum(all_edges_to_pos) <= 1,name=\
                                           "ancestor_end_edge_constraint-%s-%s"%pos_id)
                 # C - sum(edges) = position
                 self.m.addConstr(quicksum(all_edges_to_pos) == self.positions[pos_id],\
                                             name = "ancestor_edge_node_recon_constraint-%s-%s"%pos_id)


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

                    # check if it is a extant node or not
                    try:
                        node_neighbor_pos_var = int(self.extant_data[node_neighbor_item][pos])
                    except:
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

                    # difference constraints
                    diff_id  = (node,node_neighbor_item,pos)
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


    def save(self):

        # Optimize
        self.total_objective = sum([o for o in self.objective])
        self.m.setObjective(self.total_objective, GRB.MINIMIZE)
        self.m.update()

        self.m.write((self.mip_model_file))



def main():
    folder_location         = '/Users/sanjanatule/Documents/uq/Projects/MIPIndel/data/'
    #folder_location         = '/media/WorkingSpace/Share/mipindel/data/'

    ## Sample tree 1
    tree_name               =   'DHAD_1612' #'Anthony' #'DHAD_1612' #  #'CYP2U_359'#'CYP2U_165'
    mip_ancestor_fasta_file = folder_location + tree_name + '/mip_ancestor_indel.fasta'
    mip_model_file = folder_location + tree_name + '/mip_model.lp'

    # read input Dataset
    with open(folder_location + tree_name + '/ancestor_info.pkl', 'rb') as f:
        ancestor_data = pickle.load(f)

    with open(folder_location + tree_name + '/neighbor_dict.pkl', 'rb') as f:
        neighbor_dict = pickle.load(f)

    with open(folder_location + tree_name + '/extant_data.pkl', 'rb') as f:
        extant_data = pickle.load(f)


    # Run MIP
    start = time.time()
    print("Start Time:",start)

    # initialise the class
    PyTree = PhyloTreeMIP(extant_data,ancestor_data,tree_name,neighbor_dict,mip_ancestor_fasta_file,mip_model_file)

    # variable and position constraints for extants
    #print("Adding Extant Constraints")
    #PyTree.add_pos_constraints_extants()
    # variable and position constraints for ancestors
    print("Adding Ancestor Constraints")
    PyTree.add_pos_constraints_ancestors()
    # variable and edge constraints for ancestors
    print("Adding Edges Constraints")
    PyTree.add_edge_constraints_ancestors()
    # position difference constraints
    print("Adding Penalty Constraints")
    PyTree.add_penalty_constraint()

    print('Saving MIP Model')
    PyTree.save()

if __name__ == "__main__":
  main()

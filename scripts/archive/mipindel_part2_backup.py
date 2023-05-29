''' Create and run MIP model using pickle files created in part 1 '''

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


    # penalty constraint
    def add_penalty_constraint(self):
        ''' function to add penalty constraints for whole tree  '''

        # difference constraints
        for node,node_neighbor in self.neighbor_dict.items():
            for node_neighbor_item in node_neighbor:

                # V - add position difference variable
                node_pos_var = self.ancestorsequence[node]

                # check if it is a extant node or not
                if node_neighbor_item in self.extant_list:
                    node_neighbor_pos_var = [int(i) for i in self.extant_data[node_neighbor_item]]
                else:
                    node_neighbor_pos_var = self.ancestorsequence[node_neighbor_item]

                # V - create variables
                nn_pair               = (node,node_neighbor_item)
                diff_pos              = self.m.addVars(self.sequence_range,vtype=GRB.BINARY)
                pen                   = self.m.addVars(self.sequence_range,vtype=GRB.BINARY)
                self.diff[nn_pair]    = diff_pos
                self.penalty[nn_pair] = pen


                for pos in range(1,self.sequence_length - 1):
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
        self.m.Params.LogFile = "mip_model.log"

        # Optimize
        self.total_objective = quicksum(self.objective)
        self.m.setObjective(self.total_objective, GRB.MINIMIZE)
        self.m.update()

        #self.m.write(('mip_model_' + self.tree_name + '.lp'))
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

def main():
    folder_location         = '/Users/sanjanatule/Documents/uq/Projects/MIPIndel/data/'
    #folder_location         = '/media/WorkingSpace/Share/mipindel/data/'

    ## Sample tree 1
    for tree_name in ['st1']: #'CYP2U_359','CYP2U_165','DHAD_1612','ALS','anthony']:

        print("Running MIP for tree-",tree_name)

        mip_ancestor_fasta_file = folder_location + tree_name + '/mip_ancestor_indel.fasta'
        mip_model_file = folder_location + tree_name + '/mip_model.lp'

        # read input Dataset
        with open(folder_location + tree_name + '/ancestor_info.pkl', 'rb') as f:
            ancestor_data = pickle.load(f)

        with open(folder_location + tree_name + '/neighbor_dict.pkl', 'rb') as f:
            neighbor_dict = pickle.load(f)

        with open(folder_location + tree_name + '/extant_data.pkl', 'rb') as f:
            extant_data = pickle.load(f)

        # Create the MIP model
        start = time.time()
        print("Start Time:",start)
        # initialise the class
        PyTree = PhyloTreeMIP(extant_data,ancestor_data,tree_name,neighbor_dict,mip_ancestor_fasta_file,mip_model_file)
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

if __name__ == "__main__":
  main()

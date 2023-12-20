''' MIP Model build and optimize '''

# libraries
from collections import defaultdict
from ete3 import Tree
import numpy as np
from pysam import FastaFile,FastxFile
import pickle
import sys
import getopt
import gurobipy as gp
from gurobipy import abs_,quicksum,GRB
import time
import pickle

# class for mip model building and training
class PhyloTreeMIP:
    def __init__(self,extant_data,ancestor_data,neighbor_dict,alpha_p):

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
        self.alpha_p = alpha_p

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
                            
                            self.m.addConstr(pen[pos] >= node_neighbor_pos_var[pos-1] + (1 - node_pos_var[pos-1]) + (1 - node_neighbor_pos_var[pos]) + node_pos_var[pos] -3)
                            self.m.addConstr(pen[pos] >= node_neighbor_pos_var[pos] + (1 - node_pos_var[pos]) + (1 - node_neighbor_pos_var[pos-1]) + node_pos_var[pos-1] -3)

                        # O - add difference to the objective
                        self.objective.append(self.diff[(node,node_neighbor_item,pos)])

                        # O - add penalty to the objective
                        self.objective.append(self.alpha_p * pen[pos])

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
                            self.m.addConstr(pen[pos] >= node_neighbor_pos_var[pos-1] + (1 - node_pos_var[pos-1]) + (1 - node_neighbor_pos_var[pos]) + node_pos_var[pos] -3)
                            self.m.addConstr(pen[pos] >= node_neighbor_pos_var[pos] + (1 - node_pos_var[pos]) + (1 - node_neighbor_pos_var[pos-1]) + node_pos_var[pos-1] -3)


                        # O - add difference to the objective
                        self.objective.append(diff_pos[pos])

                        # O - add penalty to the objective
                        self.objective.append(self.alpha_p * pen[pos])

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
        #self.m.write('mip_log.lp')
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

    def output_fasta(self,all_node_paths,output_folder_location):
        # convert output file to FASTA file
        with open(output_folder_location + self.mip_ancestor_fasta_file,mode='w') as fout:
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

        # understand the objective score
        # score_dict = {}
        # overall_score = 0
        # for node,node_neighbor in self.neighbor_dict.items():
        #     for node_neighbor_item in node_neighbor:
        #         total_score_1 = 0
        #         total_score_2 = 0
        #         pen_id  = (node,node_neighbor_item)
        #         diff_id = (node,node_neighbor_item)
                
        #         for pos in range(1,self.sequence_length - 1): # penalty start from 1st position only
                    
        #             total_score_1 = total_score_1 + 2 * int(self.penalty[pen_id][pos].X)
        #             total_score_2 = total_score_2 + int(self.diff[diff_id][pos].X)
                
        #         print(f"Sequence for node {node} is {all_node_paths[node]}")
        #         print(f"Sequence for node {node_neighbor_item} is {all_node_paths[node_neighbor_item]}")
        #         print(f"Difference Score between {node} and {node_neighbor_item} is {total_score_2}")
        #         print(f"Penalty Score between {node} and {node_neighbor_item} is {total_score_1}")
        #         print(f"Total score is {total_score_1 + total_score_2}")

        return all_node_paths

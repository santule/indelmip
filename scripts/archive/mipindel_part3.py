''' script to test load mip model and run'''
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
    def __init__(self,lp_model,mip_ancestor_fasta_file):
        self.mip_ancestor_fasta_file = mip_ancestor_fasta_file

        # 2 - create a new model
        self.m = gp.read(lp_model)

    def train(self,n_threads,time_out):
        # Params
        self.m.Params.Threads = n_threads
        self.m.Params.TimeLimit = time_out*60
        self.m.Params.LogToConsole = 0
        self.m.Params.Degenmoves=0

        # Optimize
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


def main():

    folder_location         = '/Users/sanjanatule/Documents/uq/Projects/MIPIndel/data/'
    #folder_location         = '/media/WorkingSpace/Share/mipindel/data/'

    ## Sample tree 1
    tree_name               = 'DHAD_1612' #'Anthony' #'DHAD_1612' #'CYP2U_359' #'CYP2U_165'
    mip_ancestor_fasta_file = folder_location + tree_name + '/mip_ancestor_indel.fasta'

    # read input Dataset
    lp_model_file = folder_location + tree_name + '/mip_model.lp'


    # Run MIP
    start = time.time()
    print("Start Time:",start)
    n_threads = 1
    time_out = 60 # 1 hour

    # initialise the class
    PyTree = PhyloTreeMIP(lp_model_file,mip_ancestor_fasta_file)

    # train the model
    is_sat = PyTree.train(n_threads, time_out)
    print("is_sat",is_sat)
    total_time = ((time.time()-start))
    print("-----------------------------")
    print("Total time = %0.2f[m]"%total_time)
    info = PyTree.get_info()
    info["total_time"] = total_time
    info["is_sat"] = is_sat
    print("info",info)

if __name__ == "__main__":
  main()

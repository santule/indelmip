''' Main program for running indel inference for phylogenetic tree - input is phylogenetic tree and alignment fasta file '''
''' Output is indel fasta file '''


# libraries
from collections import defaultdict
import numpy as np
import pickle
import sys
import getopt
import time
import pickle

# import model and utils file
import model_altopt
import utils

# process input data for MIP
def data_processing(extant_sequence_file,nwk_file_path,folder_location):

    #folder_location = '/Users/sanjanatule/Documents/uq/Projects/MIPIndel/data/st1/'

    # prepare input Dataset
    print("Processing Input Files")
    MIPIndel      = utils.IndelsInfo(extant_sequence_file,nwk_file_path,folder_location) # class
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

# build mip model and train
def mip_processing(extant_data,ancestor_data,neighbor_dict,optimal_constraint):
    # Create the MIP model
    start = time.time()
    print("Start Time:",start)
    # initialise the class
    PyTree = model_altopt.PhyloTreeMIP(extant_data,ancestor_data,neighbor_dict)
    # variable position constraints for ancestors
    print("Adding Ancestor Constraints")
    PyTree.add_pos_constraints_ancestors()
    # edge constraints for ancestors
    print("Adding Edges Constraints")
    PyTree.add_edge_constraints_ancestors()
    # position difference constraints
    print("Adding Penalty Constraints")
    PyTree.add_penalty_constraint()
    # remove optimal solution
    print(f"optimal_constraint {optimal_constraint}")
    if optimal_constraint == 'Y':
        with open('optimal_edges.pkl','rb') as f:
            optimal_edge_list = pickle.load(f)
        print("Remove Optimal Solution")
        PyTree.remove_optimal_constraint(optimal_edge_list)

    total_time = ((time.time()-start)/60)
    print("-----------------------------")
    print("Total time to create model = %0.2f[mins]"%total_time)

    # Run MIP
    print('Training MIP Model')
    start = time.time()
    print("Start Time:",start)
    n_threads = 1
    time_out  = 60 # 1 hour

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
        all_node_paths,optimal_edges_list = PyTree.get_solution()
        # save the optimal edge solution in the .pkl file
        with open('optimal_edges.pkl','wb') as f:
            pickle.dump(optimal_edges_list,f)
        PyTree.output_fasta(all_node_paths)
        print(f"MIP indel output file mip_ancestor_indel.fasta created and saved")
    else:
        print("Did not find any satisfactory solution to the model")

# main function
def main(extant_sequence_file,nwk_file_path,output_folder_location,optimal_constraint):
    # 1 - process the input data for MIP
    extant_data,ancestor_data,neighbor_dict = data_processing(extant_sequence_file,nwk_file_path,output_folder_location)
    # 2 - build model
    mip_processing(extant_data,ancestor_data,neighbor_dict,optimal_constraint)

if __name__ == "__main__":
  main(extant_sequence_file,nwk_file_path,output_folder_location,optimal_constraint)

''' python script to convert indel mip output to the grasp json ASR.json'''

import json
from pysam import FastaFile,FastxFile
from ete3 import Tree
import sys

def fasta_to_json():

    indel_fasta_solution_file = 'mip_ancestor_indel2.fasta'
    asr_json_file = 'ASR.json'

    # read the asr json file
    f = open(asr_json_file)
    asr_json_data = json.load(f)

    # create edge repeat structure
    edges_repeat = {"Recip": True,"Backward": True,"Forward": True,"Weight": 0}

    # read the indel inference from mip
    ancestor_indel_info = {}
    with FastxFile(indel_fasta_solution_file) as fh:
        for entry in fh:
            ancestor_indel_info[entry.name] = entry.sequence

    # read each ancestor and change it
    for a in asr_json_data["Ancestors"]:
        # date change
        #a["GRASP_version"] = "1900-01-01"
        a["Directed"] = True
        a["Terminated"] = True
        ancestor_name  = "N" + a["Name"]
        ancestor_indel = ancestor_indel_info[ancestor_name]

        # read the mip indel and change the values
        # change the Indices
        indices = [e_idx - 1 for e_idx,e_val in enumerate(list(ancestor_indel)) if e_val == '1']
        a["Indices"] = indices[1:-1] # do not include first and last node as mip has start/end node
        # changes the number of nodes to be same as the indices
        a["Nodes"] = [{} for i in range(len(a["Indices"]))]

        # change the Edgeindices and adjacent structure
        last_node = a["Indices"][-1]
        all_edges = []
        adjacent  = []
        for i_dx,i_val in enumerate(indices[:-1]):
            all_edges.append([i_val, indices[i_dx + 1]])
            if i_val != -1: # starting node
                if i_val != last_node: # last node
                    adjacent.append([indices[i_dx + 1]])
                else:
                    adjacent.append([])

        a["Edgeindices"] = all_edges
        a["Adjacent"]    = adjacent

        # edges for all edges in the edge indices
        a["Edges"] = [edges_repeat for i in range(len(a["Edgeindices"]))]

    # Dump the output in the file
    with open('ASR_MIP2.json', 'w') as outfile:
        json.dump(asr_json_data, outfile)


if __name__ == "__main__":
  fasta_to_json()

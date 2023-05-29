from pysam import FastaFile,FastxFile
import sys
import json

def read_sequences(ancestor_file,extant_file,output_json):

    #output_json = '/Users/sanjanatule/Documents/uq/Projects/PreferredPath/data/tree_annotations.json'
    #sequence_dict = {}
    json_dict = {"treename": "Demo","nodes":{},"important":["root"]}


    # read ancestor file
    with FastxFile(ancestor_file) as fh:
        for entry in fh:
            #sequence_dict[entry.name] = str(entry.sequence)
            json_dict["nodes"][entry.name] = {}
            json_dict["nodes"][entry.name]["SEQ"] = str(entry.sequence)
            #str(entry.sequence)

    # read extant file
    with FastxFile(extant_file) as fh:
        for entry in fh:
            #sequence_dict[entry.name] = str(entry.sequence)
            json_dict["nodes"][entry.name] = {}
            json_dict["nodes"][entry.name]["SEQ"] = str(entry.sequence)

    # make special entry for root as root is called N0 in grasp
    json_dict["nodes"]['root'] = json_dict["nodes"]['N0']


    

    # create json file
    with open(output_json, mode='w') as fout:
        json.dump(json_dict,fout)
        fout.close()


def main():
    # input file and output file
    ancestor_file = sys.argv[1]
    extant_file   = sys.argv[2]
    output_json   = sys.argv[3]
    read_sequences(ancestor_file,extant_file,output_json)


if __name__ == "__main__":
    main()

''' Python script to convert GRASP output into binary .fasta file '''

import json
from pysam import FastaFile,FastxFile
import sys
import getopt

''' function to convert GRASP json output to binary fasta file '''
def convert_grasp_json_fasta_file(**kwargs):

    # open and parse json file
    f = open(kwargs['json_file_name'])
    json_data = json.load(f)

    # process json
    ancestor_indel_list = []
    ancestor_name = []

    for a_info in json_data['Ancestors']: # parse ancestor
        an_name = a_info['Name']
        an_indel_string = ['0'] * a_info['Size'] # default string with size
        for ei in a_info['Edgeindices']: # edge indices process for an ancestor
            from_pos,to_pos = ei[0],ei[1]
            if from_pos in a_info['Indices']: # start and end node needs to be ignored
                an_indel_string[from_pos] = '1'
            if to_pos in a_info['Indices']:
                an_indel_string[to_pos] = '1'
        an_indel_string = ''.join(an_indel_string)
        ancestor_indel_list.append(an_indel_string)
        ancestor_name.append(an_name)

    # save the solution in the .fasta file
    with open(kwargs['output_folder_location'] + '/' + kwargs['grasp_indel_method'] + '_grasp_all_indel.fasta',mode = 'w') as fw:
        for name,indel in zip(ancestor_name,ancestor_indel_list):
            fw.write('>' + str(name) + '\n')
            fw.write(str(indel) + '\n')
    print(f"Successfully created the indel fasta file - {kwargs['grasp_indel_method'] + '_grasp_all_indel.fasta'}")


''' function to convert GRASP alphabet fasta output to binary fasta file '''
def convert_grasp_fasta_to_binary(**kwargs):

    # process json
    ancestor_indel_list = []
    ancestor_name = []

    with FastxFile(kwargs['output_folder_location'] + '/' + kwargs['fasta_file_name']) as fh:
        for entry in fh:
            # add start and end string to the sequence
            seq_name = entry.name

            # binarise sequences
            seq_binary   = [0 if s == '-' else 1 for s in entry.sequence]
            seq_binary   = ''.join([str(x) for x in seq_binary])

            # add to the list
            ancestor_name.append(seq_name)
            ancestor_indel_list.append(seq_binary)

    # save the solution in the .fasta file
    with open(kwargs['output_folder_location'] + '/' + kwargs['grasp_indel_method'] + '_grasp_all_indel.fasta',mode = 'w') as fw:
        for name,indel in zip(ancestor_name,ancestor_indel_list):
            fw.write('>' + str(name) + '\n')
            fw.write(str(indel) + '\n')
    print(f"Successfully created the indel fasta file - {kwargs['grasp_indel_method'] + '_grasp_all_indel.fasta'}")

''' function to combine ancestor and extant .fasta or .aln files '''
def combine_extant_ancestor_data(**kwargs):

    extant_info = ancestor_info = ''
    with open(kwargs['extant_fasta_file']) as file1:
        extant_info = file1.read()

    with open(kwargs['ancestor_fasta_file']) as file2:
        ancestor_info = file2.read()

    extant_info += '\n'
    extant_info += ancestor_info

    with open(kwargs['output_folder_location'] + '/' + kwargs['grasp_indel_method'] + '_grasp_all.fasta' , 'w') as file:
        file.write(extant_info)
    print(f"Successfully created the fasta file - {kwargs['grasp_indel_method'] + '_grasp_all.fasta'}")


# main function
def main(ancestor_fasta_file_name,extant_fasta_file_name,grasp_indel_method,data_folder):

        pr = extant_fasta_file_name.split("/")[0]
        print(pr)
        kwargs_parameter = {'output_folder_location':data_folder,'extant_fasta_file':extant_fasta_file_name,\
                            'ancestor_fasta_file':ancestor_fasta_file_name ,'grasp_indel_method':grasp_indel_method}
        combine_extant_ancestor_data(**kwargs_parameter)
        combined_file = grasp_indel_method + '_grasp_all.fasta'
        kwargs_parameter = {'output_folder_location':data_folder,'fasta_file_name':combined_file,\
                            'grasp_indel_method':grasp_indel_method}
        convert_grasp_fasta_to_binary(**kwargs_parameter)


if __name__ == "__main__":
  main(ancestor_fasta_file_name,extant_fasta_file_name,grasp_indel_method,data_folder)

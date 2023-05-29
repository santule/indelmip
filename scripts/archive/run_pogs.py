import json
import glob
import numpy as np
import sys
import csv
import os
import re


# function to count the number of paths in a POG
def count_path(a):
    a = a + a.T    #add up the transpose
    a = np.clip(a,0,1)
    a = np.triu(a) #only the upper triangle

    nodes = a.shape[0]
    dp = [0] * nodes
    dp[nodes - 1]= 1 #last node


    for i in range(nodes - 1, -1, -1):

      # find all neighbours of the current node. ( i.e. nodes going out from current node)
      # index in the array where the current row has 1
        neighbour_nodes = np.where (a[i] == 1)[0]

        for j in neighbour_nodes:
            dp[i] = dp[i] + dp[j]

    return(dp[0])

def count_total_pogs(folder):

    # read all .dot files from the directory
    total_sequences_list = []
    all_exp_results = []

    path = folder + '/pogs.json'
    instance_name = path.split('/')[-2]

    # read the json file
    with open(path, 'r') as j:
        pog_all_data = json.loads(j.read())

        # read all ancestors
        for pog_data in pog_all_data['Ancestors']:
            node_name = 'N' + pog_data['Name']

            #print("node_name",node_name)

            # read that node's data
            nodes = pog_data['Size'] + 2

            # create numpy zero matrix
            mat = np.zeros(shape=(nodes,nodes))

            # Edges from special Start node to the start nodes
            for s in pog_data['Starts']:
                mat[0,s + 1] = 1

            # Edges from last node to the special End node
            for e in pog_data['Ends']:
                mat[e + 1,nodes-1] = 1

            # create the adjency matrix for all nodes except from special node start
            for e in pog_data['Edges']:
                #print("e[From]={} and e[To]={}".format(e['From'],e['To']))
                mat[e['From'] + 1,e['To'] + 1] = 1

            #print("mat",mat)
            # count the number of viable paths from the POG
            total_sequences = count_path(mat)
            #print("total_sequences",total_sequences)

            # append path to the .csv file
            total_sequences_list.append([instance_name,node_name,total_sequences])

    # all_exp_results.append((instance_name , total_sequences_list))
    return total_sequences_list


def main():
    # input file and output file
    input_sub_folder = sys.argv[1]
    all_exp_results = count_total_pogs(input_sub_folder)

    # save to .csv files
    with open(input_sub_folder + '/pogs_count.csv', 'w',encoding='UTF8', newline='') as f:

        # using csv.writer method from CSV package
        write = csv.writer(f)
        write.writerows(all_exp_results)
        f.close()


if __name__ == "__main__":
    main()

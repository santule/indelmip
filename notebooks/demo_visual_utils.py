import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os, glob
from matplotlib.lines import Line2D 
from pysam import FastaFile,FastxFile
from ete3 import Tree
import numpy as np

def get_indel_events(str1,str2):
    dis = 0
    prev_dis = 0

    for i in range(0,len(str1)):
        curr_dis = int(str1[i]) - int(str2[i])
        
        if curr_dis != 0 and curr_dis != prev_dis:
            dis += 1
        prev_dis = curr_dis
        
    return dis
def get_whole_tree_indel_events(method_fasta_file,nwk_file_path,col,method_name):
    level = 0
    extant_list = []
    
    # sequence info
    sequences_fasta_info = FastaFile(method_fasta_file)
    
    # tree file
    tree_file = open(nwk_file_path,"r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)
    
    for n in tree.traverse():
        if n.up is not None: # root node
            n.add_features(level = n.up.level + 1)
            seq_name = n.name
            indel_sequence_curr_level = sequences_fasta_info.fetch(n.name)
            indel_sequence_up_level   = sequences_fasta_info.fetch(n.up.name)
            total_indel_events = get_indel_events(indel_sequence_up_level,indel_sequence_curr_level)
            n.add_features(indel_events = total_indel_events + n.up.indel_events)
        else:
            n.add_features(level = level)
            n.add_features(indel_events = 0)

        if n.is_leaf() == True:
            extant_list.append(n.name)

    tree_exp_indel_event = []
    for ext in extant_list:
        x_node_trav = []
        y_node_trav = []
        node = tree.search_nodes(name=ext)[0]
        while node:
            x_node_trav.append(node.level)
            y_node_trav.append(node.indel_events)
            node = node.up
        tree_exp_indel_event.append([x_node_trav,y_node_trav])
        
    return tree_exp_indel_event

def combined_tree_visual(path_link):

    nwk_file_path           = path_link + '/psp_ancestors.nwk'
    mip_indel_file          = path_link + '/mip_ancestor_indel.fasta'  # indels from MIP
    grasp_indel_file        = path_link + '/bep_grasp_all_indel.fasta' # indels from grasp
    sicp_indel_file         = path_link + '/sicp_grasp_all_indel.fasta' # indels from grasp
    psp_indel_file          = path_link + '/psp_grasp_all_indel.fasta' # indels from grasp

    bep_indel_event  = get_whole_tree_indel_events(grasp_indel_file,nwk_file_path,'#0F52BA',"BEP")
    psp_indel_event  = get_whole_tree_indel_events(psp_indel_file,  nwk_file_path, '#DD3300',"PSP")
    sicp_indel_event = get_whole_tree_indel_events(sicp_indel_file, nwk_file_path, '#00DD03',"SICP")
    mip_indel_event  = get_whole_tree_indel_events(mip_indel_file,  nwk_file_path, '#58181F',"MIP")

    plt.rcParams["figure.figsize"] = [5,5] 
    plt.rcParams["figure.dpi"] = 1000
    plt.ylim(0, 140)
    plt.xlim(0,30)
    plt.xlabel("Tree Level")
    plt.ylabel("Indel Events")


    legend_elements = [Line2D([0], [0], color='#0F52BA', lw=2, label='BEP'),
                    Line2D([0], [0], color='#DD3300', lw=2, label='PSP'),
                    Line2D([0], [0], color='#00DD03', lw=2, label='SICP'),
                    Line2D([0], [0], color='#58181F', lw=2, label='MIP')]

    for pts in bep_indel_event:
        plt.plot(pts[0],pts[1], '-o',color='#0F52BA',alpha=0.3,linewidth=0.5,markersize=1)
    for pts in psp_indel_event:
        plt.plot(pts[0],pts[1], '-o',color='#DD3300',alpha=0.3,linewidth=0.5,markersize=1)
    for pts in sicp_indel_event:
        plt.plot(pts[0],pts[1], '-o',color='#00DD03',alpha=0.3,linewidth=0.5,markersize=1)
    for pts in mip_indel_event:
        plt.plot(pts[0],pts[1], '-o',color='#58181F',alpha=0.6,linewidth=0.5,markersize=1)

    plt.legend(handles=legend_elements, loc='upper left')
    plt.show()

def eval_metrics_plot(path_link):
    temp_list = []
    # indel score
    os.chdir(path_link)
    for fname in glob.glob("*indscore*"):
        with open(fname) as infile:
            for line in infile:
                temp_list.append(line.split(','))
    df_is = pd.DataFrame(temp_list,columns=['Method','Indel score'])
    
    # out of distibution
    temp_list = []
    # indel score
    os.chdir(path_link)
    for fname in glob.glob("*out_dist_percent*"):
        with open(fname) as infile:
            for line in infile:
                temp_list.append(line.split(','))
    df_od = pd.DataFrame(temp_list,columns=['Method','Percentage of ancestors with out of distribution pattern'])
    
    # 3 mutations
    temp_list = []
    # indel score
    os.chdir(path_link)
    for fname in glob.glob("*ancestors_with_3_mut*"):
        with open(fname) as infile:
            for line in infile:
                temp_list.append(line.split(','))
    df_3m = pd.DataFrame(temp_list,columns=['Method','Percentage of ancestors with three mutations away'])
    
    # plot them
    plt.figure(num=None, figsize=(16, 6))
    font = {'family' : 'Times New Roman',
            'size'   : 10}

    plt.rc('font', **font)

    plt.subplot(1,3,1)
    plt.bar(df_is['Method'],df_is['Indel score'].astype(int))
    plt.ylabel("Indel score")
    
    plt.subplot(1,3,2)
    plt.bar(df_od['Method'],df_od['Percentage of ancestors with out of distribution pattern'].astype(float))
    plt.ylabel("Percentage of ancestors with out of distribution pattern")
    
    plt.subplot(1,3,3)
    plt.bar(df_3m['Method'],df_3m['Percentage of ancestors with three mutations away'].astype(float))
    plt.ylabel("Percentage of ancestors with three mutations away")
    
    plt.show()
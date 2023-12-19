'''
Python scripts to run mip inference method for different hyperparameter alpha values
'''

import os
import subprocess
import sys
import getopt
import score_objective,score_indel,score_distribution,score_ancestor_3_mutation_away

# main function
def exp_hyperparameter(dataset_folder):
    exp_results = []
    for alpha in [0,1,2,3,4,5,6,7,8]:
        print(f"Processing dataset {dataset_folder} with alpha parameter {alpha}")

        align_file = dataset_folder + 'extants.aln'
        output_folder = dataset_folder

        print(f"1 - Running MIP")
        tree_file = dataset_folder + 'psp_ancestors.nwk'
        log_file  = dataset_folder + 'mip_inference_log.txt'
        cmd = "python run_mipindel.py -a {align_file}\
            -n {tree_file} -p {alpha_parameter} -o {output_folder}\
            > {log_file}".format(align_file=align_file,\
            tree_file=tree_file,log_file=log_file,alpha_parameter = alpha,output_folder=dataset_folder)
        
        print(cmd)
        subprocess.run(cmd,shell=True)

        ## evaluation
         # Objective Scoring
        print("2 - RUNNING OBJECTIVE SCORING")
        objscore = score_objective.score_fn(tree_file,dataset_folder + 'mip_ancestor_indel.fasta')

        # Indel Scoring
        print("3 - RUNNING INDEL SCORING")
        indscore = score_indel.score_fn(tree_file,dataset_folder + 'mip_ancestor_indel.fasta')

        exp_results.append([alpha,objscore,indscore])

    return exp_results
       

if __name__ == "__main__":
    dataset_folder    = sys.argv[1]   
    exp_hyperparameter(dataset_folder)

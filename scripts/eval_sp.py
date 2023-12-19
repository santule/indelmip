'''
Python scripts to compare Srinr and Practer to the MIP Indel inference method.
'''

import os
import subprocess
import sys
from os import walk
import time


# main function
def main(data_folder):

    # Run evaluation scripts
    for (sub_folder, _, _) in walk(data_folder):
        if sub_folder != data_folder:
            sub_folder = sub_folder + '/'

            print(f"******** Processing protein family {sub_folder}")
            align_file = sub_folder +  'extants.aln'
            tree_file  = sub_folder + 'input_tree.nwk'
            output_folder = sub_folder

            ## grasp inference to get the ancestral nodes names
            print(f"Processing PSP Inference")
            log_file   = sub_folder + 'psp_inference_log.txt'
            cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method PSP --onlyindel --save-all -pre 'psp' --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)

            ## MIP Inference
            print(f"Running MIP")
            tree_file = sub_folder + 'psp_ancestors.nwk'
            log_file = sub_folder + 'mip_inference_log.txt'
            alpha = 2
            start_time = time.time()
            print(f"Time start :{start_time}")
            cmd = "python run_mipindel.py -a {align_file} -n {tree_file} -p {alpha_parameter} -o {output_folder} > {log_file}"\
                .format(align_file=align_file,tree_file=tree_file,log_file=log_file,alpha_parameter = alpha,output_folder=output_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)
            print(f"Time finish :{time.time()}")
            print(f"Total time in secs: {((time.time()-start_time))}")
            print("Finished running.")

            ## MIP evaluation
            tree_file = sub_folder + 'psp_ancestors.nwk'
            eval_log_file = sub_folder + 'mip_evaluate_log.txt'
            cmd = "python metrics_mip.py  -f {pr_folder} -e {align_file} -n {tree_file} > {eval_log_file}".format\
                (align_file = align_file,tree_file=tree_file,\
                 eval_log_file=eval_log_file,pr_folder = sub_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)

            # SP inference
            print("Running Sr and Pt algorithm")
            start_time = time.time()
            print(f"Time start :{start_time}")
            cmd = "python run_spindel.py -e {align_file} -n {tree_file}".format(align_file = align_file,tree_file=tree_file)
            print(cmd)
            subprocess.run(cmd,shell=True)
            print(f"Time finish :{time.time()}")
            print(f"Total time in secs: {((time.time()-start_time))}")
            print("Finished running.")

            # SP evaluation
            print("Running evaluation")
            eval_log_file = sub_folder + 'sp_evaluate_log.txt'
            cmd = "python metrics_sp.py -f {pr_folder} -e {align_file} -n {tree_file} > {eval_log_file}"\
                .format(align_file = align_file,\
                        tree_file=tree_file,eval_log_file=eval_log_file,pr_folder = sub_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)
            print(f"************************\n")


if __name__ == "__main__":
    data_folder = sys.argv[1]
    main(data_folder)
    print("Done experiment.")

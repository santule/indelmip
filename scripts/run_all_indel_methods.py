'''
Python scripts to run all indel inference methods. 
For MIP, branchpoints needs to be named, hence psp is run first.
'''

import os
import subprocess
import sys
import getopt

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python run_all_indel_methods.py -a alignment file -n newick tree file -o output file location -b y -m y -p y -s y')
    exit()

# main function
def eval_indel():

    align_file = None
    tree_file = None
    output_folder  = None
    run_mip   = None
    run_bep   = None
    run_sicp  = None
    run_psp   = None

    argv    = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:n:o:b:m:p:s:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-a']:
            align_file = arg
        elif opt in ['-n']:
            tree_file = arg
        elif opt in ['-o']:
            output_folder = arg
        elif opt in ['-b']:
            run_bep  = arg
        elif opt in ['-m']:
            run_mip  = arg
        elif opt in ['-p']:
            run_psp  = arg
        elif opt in ['-s']:
            run_sicp = arg

    if align_file is None or tree_file is None or output_folder is None or run_bep is None or run_mip is None or run_psp is None or run_sicp is None:
        help()

    print(f"Processing Indel Inference") 
    # Run evaluation scripts
    if run_bep == 'y':
        # BEP
        ## Indel Inference
        print(f"Processing BEP Inference")
        log_file   = output_folder + 'bep_inference_log.txt'
        cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method BEP --onlyindel --save-all --threads 1 --orphans -pre 'bep' -o {output_folder}> {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

        ## AA Inference
        log_file   = output_folder +  'bep_inference_full_log.txt'
        cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method BEP --orphans -pre 'bep' -o {output_folder} > {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

        ## Evaluation
        grasp_tree = output_folder + 'bep_ancestors.nwk'
        eval_log_file = output_folder + 'bep_evaluate_log.txt'
        grasp_ancestor_file = output_folder + 'bep_ancestors.fa'
        cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"bep\" > {eval_log_file}".\
            format(grasp_ancestor_file=grasp_ancestor_file,\
                    align_file = align_file,grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder = output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

    if run_psp == 'y':
        # PSP
        ## Inference
        print(f"Processing PSP Inference")
        log_file   = output_folder +  'psp_inference_log.txt'
        cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method PSP --onlyindel --save-all -pre 'psp' --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

        log_file   = output_folder +  'psp_inference_full_log.txt'
        cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method PSP -pre 'psp' --orphans -o {output_folder} > {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

        ## evaluation
        grasp_tree = output_folder  + 'psp_ancestors.nwk'
        eval_log_file = output_folder +  'psp_evaluate_log.txt'
        grasp_ancestor_file = output_folder + 'psp_ancestors.fa'
        cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"psp\" > {eval_log_file}"\
            .format(grasp_ancestor_file=grasp_ancestor_file,align_file = align_file,\
                    grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder = output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

    if run_sicp == 'y':
        # SICP
        ## Inference
        print(f"Processing SICP Inference")
        log_file   = output_folder + 'sicp_inference_log.txt'
        cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method SICP --onlyindel --save-all -pre 'sicp'  --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

        log_file   = output_folder + 'sicp_inference_full_log.txt'
        cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method SICP -pre 'sicp' --orphans -o {output_folder} > {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

        ## evaluation
        grasp_tree = output_folder + 'sicp_ancestors.nwk'
        eval_log_file = output_folder + 'sicp_evaluate_log.txt'
        grasp_ancestor_file = output_folder + 'sicp_ancestors.fa'
        cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"sicp\" > {eval_log_file}".\
            format(grasp_ancestor_file=grasp_ancestor_file,\
                    align_file = align_file,grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder =output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

    if run_mip == 'y':
        # MIP
        ## Inference
        print(f"Processing MIP Inference")
        tree_file = output_folder + 'psp_ancestors.nwk'
        log_file  = output_folder + 'mip_inference_log.txt'
        alpha = 2
        cmd = "python run_mipindel.py -a {align_file} -n {tree_file} -p {alpha_parameter} -o {output_folder} >\
            {log_file}".format(align_file=align_file,tree_file=tree_file,log_file=log_file,alpha_parameter = alpha,output_folder=output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)

        ## evaluation
        tree_file = output_folder + 'psp_ancestors.nwk'
        eval_log_file = output_folder + 'mip_evaluate_log.txt'
        cmd = "python metrics_mip.py -f {pr_folder} -a {align_file} -n {tree_file} > {eval_log_file}".\
            format(align_file = align_file,tree_file=tree_file,eval_log_file=eval_log_file,pr_folder =output_folder)
        print(cmd)
        subprocess.run(cmd,shell=True)


if __name__ == "__main__":
    eval_indel()
    print("Completed Indel Inference")

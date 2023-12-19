'''
Python scripts to run all indel inference methods. Run this script from the data folder.
The data folders contains sub folders, each for one protein family.
Sub folder expected to contain .aln and .nwk file.
For MIP, branchpoints needs to be named, hence psp is run first.
'''

import os
import subprocess
import sys
import getopt

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python run_indel_real_datasets.py -f data_folder -b y -m y -p y -s y')
    exit()

# main function
def eval_real():
    run_mip   = None
    run_bep   = None
    run_sicp  = None
    run_psp   = None
    data_folder = None

    argv    = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "f:b:m:p:s:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-b']:
            run_bep  = arg
        elif opt in ['-m']:
            run_mip  = arg
        elif opt in ['-p']:
            run_psp  = arg
        elif opt in ['-s']:
            run_sicp = arg
        elif opt in ['-f']:
            data_folder = arg

    if run_bep is None or run_mip is None or run_psp is None or run_sicp is None or data_folder is None:
        help()

    # Run evaluation scripts
    pr_families = ['CYP2U_165','B3_225','RNaseZ_243','CYP2U_359','GDH-GOx_399','DHAD_585','CYP2U_595','KARI_716','KARI_1176','ALPHA_1263','DHAD_1612','DHAD_1658','ALS_1990','RNaseZ_624','CYP_1656']
    for pr in pr_families:
        print(f"Processing protein family {pr}")
        align_file    = data_folder + pr + '/' + pr + '.aln'
        tree_file     = data_folder + pr + '/' + pr + '.nwk'
        output_folder = data_folder + pr + '/'

        if run_bep == 'y':
            if pr not in ['MBL_624']:
                # BEP
                ## Inference
                print(f"Processing BEP Inference")
                log_file   = data_folder + pr + '/' + 'bep_inference_log.txt'
                cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method BEP --onlyindel --save-all --threads 1 --orphans -pre 'bep' -o {output_folder}> {log_file}".format(\
                            align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                log_file   = data_folder + pr + '/' +  'bep_inference_full_log.txt'
                cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method BEP --orphans -pre 'bep' -o {output_folder} > {log_file}".format(\
                            align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                ## evaluation
                grasp_tree = data_folder + pr + '/' + 'bep_ancestors.nwk'
                eval_log_file = data_folder + pr + '/' + 'bep_evaluate_log.txt'
                grasp_ancestor_file = data_folder + pr + '/' + 'bep_ancestors.fa'
                cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"bep\" > {eval_log_file}".\
                    format(grasp_ancestor_file=grasp_ancestor_file,\
                           align_file = align_file,grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder = data_folder + pr + '/')
                print(cmd)
                subprocess.run(cmd,shell=True)

        if run_psp == 'y':
            # PSP
            ## Inference
            print(f"Processing PSP Inference")
            log_file   = data_folder + pr + '/' + 'psp_inference_log.txt'
            cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method PSP --onlyindel --save-all -pre 'psp' --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)

            log_file   = data_folder + pr + '/' + 'psp_inference_full_log.txt'
            cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method PSP -pre 'psp' --orphans -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)

            ## evaluation
            grasp_tree = data_folder + pr + '/' + 'psp_ancestors.nwk'
            eval_log_file = data_folder + pr + '/' + 'psp_evaluate_log.txt'
            grasp_ancestor_file = data_folder + pr + '/' + 'psp_ancestors.fa'
            cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"psp\" > {eval_log_file}"\
                .format(grasp_ancestor_file=grasp_ancestor_file,align_file = align_file,\
                        grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder = data_folder + pr + '/')
            print(cmd)
            subprocess.run(cmd,shell=True)

        if run_sicp == 'y':
            # SICP
            ## Inference
            print(f"Processing SICP Inference")
            log_file   = data_folder + pr + '/' + 'sicp_inference_log.txt'
            cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method SICP --onlyindel --save-all -pre 'sicp'  --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)

            log_file   = data_folder + pr + '/' + 'sicp_inference_full_log.txt'
            cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method SICP -pre 'sicp' --orphans -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)

            ## evaluation
            grasp_tree = data_folder + pr + '/' + 'sicp_ancestors.nwk'
            eval_log_file = data_folder + pr + '/' + 'sicp_evaluate_log.txt'
            grasp_ancestor_file = data_folder + pr + '/' + 'sicp_ancestors.fa'
            cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"sicp\" > {eval_log_file}".\
                format(grasp_ancestor_file=grasp_ancestor_file,\
                       align_file = align_file,grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder = data_folder + pr + '/')
            print(cmd)
            subprocess.run(cmd,shell=True)

        if run_mip == 'y':
            # MIP
            ## Inference
            print(f"Running MIP")
            tree_file = data_folder + pr + '/' + 'psp_ancestors.nwk'
            log_file  = data_folder + pr + '/' + 'mip_inference_log.txt'
            alpha = 2
            cmd = "python run_mipindel.py -a {align_file} -n {tree_file} -p {alpha_parameter} -o {output_folder} >\
                {log_file}".format(align_file=align_file,tree_file=tree_file,log_file=log_file,alpha_parameter = alpha,output_folder=output_folder)
            print(cmd)
            subprocess.run(cmd,shell=True)

            ## evaluation
            tree_file = data_folder + pr + '/' + 'psp_ancestors.nwk'
            eval_log_file = data_folder + pr + '/' + 'mip_evaluate_log.txt'
            cmd = "python metrics_mip.py -f {pr_folder} -e {align_file} -n {tree_file} > {eval_log_file}".\
                format(align_file = align_file,tree_file=tree_file,eval_log_file=eval_log_file,pr_folder = data_folder + pr + '/')
            print(cmd)
            subprocess.run(cmd,shell=True)


if __name__ == "__main__":
    eval_real()
    print("COMPLETED EXPERIMENTS")

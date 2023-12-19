import os
import subprocess
import sys
import getopt
from os import walk

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python run_indel_syn_datasets.py -f data_folder -b y -m y -p y -s y')
    exit()

# main function
def eval_syn():
    run_mip   = None
    run_bep   = None
    run_sicp  = None
    run_psp   = None
    data_folder = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "f:b:m:p:s:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-b']:
            run_bep = arg
        elif opt in ['-m']:
            run_mip = arg
        elif opt in ['-p']:
            run_psp = arg
        elif opt in ['-s']:
            run_sicp = arg
        elif opt in ['-f']:
            data_folder = arg

    if run_bep is None or run_mip is None or run_psp is None or run_sicp is None or data_folder is None:
        help()
    
    # run the script from travis folder...
    for (sub_folder, _, _) in walk(data_folder):
        if sub_folder != data_folder:
            sub_folder = sub_folder + '/'

            print(f"Processing synthetic protein family {sub_folder}")
            align_file = sub_folder + 'extants.aln'
            tree_file  = sub_folder + 'input_tree.nwk'
            output_folder = sub_folder 

            if run_bep == 'y':
                # BEP
                ## Inference
                print(f"Processing BEP Inference")
                log_file   = sub_folder + 'bep_inference_log.txt'
                cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method BEP --onlyindel --save-all --threads 1 --orphans -pre 'bep' -o {output_folder}> {log_file}".format(\
                            align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                log_file   = sub_folder +  'bep_inference_full_log.txt'
                cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method BEP --orphans -pre 'bep' -o {output_folder} > {log_file}".format(\
                            align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                ## evaluation
                grasp_tree = sub_folder + 'bep_ancestors.nwk'
                eval_log_file = sub_folder + 'bep_evaluate_log.txt'
                grasp_ancestor_file = sub_folder + 'bep_ancestors.fa'
                cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"bep\" >\
                    {eval_log_file}".\
                    format(grasp_ancestor_file=grasp_ancestor_file,align_file = align_file,\
                           grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder = sub_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

            if run_psp == 'y':
                # PSP
                ## Inference
                print(f"Processing PSP Inference")
                log_file   = sub_folder + 'psp_inference_log.txt'
                cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method PSP --onlyindel --save-all -pre 'psp' --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                log_file   = sub_folder + 'psp_inference_full_log.txt'
                cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method PSP -pre 'psp' --orphans -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                ## evaluation
                grasp_tree = sub_folder + 'psp_ancestors.nwk'
                eval_log_file = sub_folder + 'psp_evaluate_log.txt'
                grasp_ancestor_file = sub_folder + 'psp_ancestors.fa'
                cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"psp\" >\
                      {eval_log_file}".\
                        format(grasp_ancestor_file=grasp_ancestor_file,align_file = align_file,\
                               grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder = sub_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

            if run_sicp == 'y':
                # SICP
                ## Inference
                print(f"Processing SICP Inference")
                log_file   = sub_folder + 'sicp_inference_log.txt'
                cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method SICP --onlyindel --save-all -pre 'sicp'  --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                log_file   = sub_folder + 'sicp_inference_full_log.txt'
                cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method SICP -pre 'sicp' --orphans -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                ## evaluation
                grasp_tree = sub_folder + 'sicp_ancestors.nwk'
                eval_log_file = sub_folder + 'sicp_evaluate_log.txt'
                grasp_ancestor_file = sub_folder + 'sicp_ancestors.fa'
                cmd = "python metrics_grasp.py -f {pr_folder} -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"sicp\" >\
                      {eval_log_file}".format(grasp_ancestor_file=grasp_ancestor_file,\
                      align_file = align_file,grasp_tree=grasp_tree,eval_log_file=eval_log_file,pr_folder = sub_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

            if run_mip == 'y':
                # MIP
                ## Inference
                print(f"Running MIP")
                tree_file = sub_folder + 'psp_ancestors.nwk'
                log_file = sub_folder + 'mip_inference_log.txt'
                alpha = 2
                cmd = "python run_mipindel.py -a {align_file} -n {tree_file} -p {alpha_parameter} -o {output_folder} > {log_file}".\
                    format(align_file=align_file,tree_file=tree_file,log_file=log_file\
                           ,alpha_parameter = alpha,output_folder= output_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)

                ## evaluation
                tree_file = sub_folder + 'psp_ancestors.nwk'
                eval_log_file = sub_folder + 'mip_evaluate_log.txt'
                cmd = "python metrics_mip.py  -f {pr_folder} -e {align_file} -n {tree_file} > {eval_log_file}".\
                    format(align_file = align_file,tree_file=tree_file,\
                           eval_log_file=eval_log_file,pr_folder = sub_folder)
                print(cmd)
                subprocess.run(cmd,shell=True)


if __name__ == "__main__":
    eval_syn()
    print("COMPLETED EXPERIMENTS")

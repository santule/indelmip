import os,sys
import subprocess


# main function
def main(data_folder):
    
    # Run evaluation scripts
    for pr in ['RHYS_1263']: #ALPHA_1263','CYP_1656']:
 #   for pr in ['CYP2U_165','MBL_243','CYP2U_359','GDH-GOx_399','DHAD_585','CYP2U_595','KARI_716','KARI_1176','DHAD_1612','DHAD_1658','ALS_1990','MBL_624']: 

        print(f"Processing protein family {pr}")
        align_file = pr + '.aln'
        tree_file  = pr + '.nwk'
        output_folder = '.'

        print("Change directory")
        os.chdir(data_folder + pr)
        print(os.getcwd())

        # MIP - 1 run
        ## Inference
        print(f"Running MIP 1")
        tree_file = 'psp_ancestors.nwk'
        log_file  = 'mip_inference_log_1.txt'
        cmd = "python run_mip_altopt.py -a {align_file} -n {tree_file} -o '.' -m {opt_option} > {log_file}"\
                    .format(align_file=align_file,tree_file=tree_file,log_file=log_file,opt_option='N')
        print(cmd)
        subprocess.run(cmd,shell=True)

        ## evaluation
        tree_file     = 'psp_ancestors.nwk'
        eval_log_file = 'mip_evaluate_log_1.txt'
        cmd = "python metrics_mip_altopt1.py -e {align_file} -n {tree_file} > {eval_log_file}".format(align_file = align_file,tree_file=tree_file,eval_log_file=eval_log_file)
        print(cmd)
        subprocess.run(cmd,shell=True)

        # MIP - 2 run
        print(f"Running MIP 2")
        tree_file = 'psp_ancestors.nwk'
        log_file  = 'mip_inference_log_2.txt'
        cmd = "python run_mip_altopt.py -a {align_file} -n {tree_file} -o '.' -m {opt_option} > {log_file}"\
                .format(align_file=align_file,tree_file=tree_file,log_file=log_file,opt_option='Y')
        print(cmd)
        subprocess.run(cmd,shell=True)

        ## evaluation
        tree_file     = 'psp_ancestors.nwk'
        eval_log_file = 'mip_evaluate_log_2.txt'
        cmd = "python metrics_mip_altopt2.py -e {align_file} -n {tree_file} > {eval_log_file}".format(align_file = align_file,tree_file=tree_file,eval_log_file=eval_log_file)
        print(cmd)
        subprocess.run(cmd,shell=True)


        os.chdir('../')
        print(os.getcwd())


if __name__ == "__main__":
    data_folder = sys.argv[1]
    main(data_folder)

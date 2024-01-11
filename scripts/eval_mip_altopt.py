import os,sys,shutil
import subprocess
import getopt

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python eval_mip_altopt.py -a alignment file -n newick tree file -o output file location')
    exit()

# main function
def eval_alt_indel():

    align_file = None
    tree_file = None
    output_folder  = None

    argv    = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:n:o:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-a']:
            align_file = arg
        elif opt in ['-n']:
            tree_file = arg
        elif opt in ['-o']:
            output_folder = arg

    if align_file is None or tree_file is None or output_folder is None:
        help()

    data_folder = output_folder

    # MIP - 1 run
    ## Inference
    print(f"Running MIP - FIRST INDEL SOLUTION")
    log_file  = output_folder + 'mip_inference_log_1.txt'
    cmd = "python run_mipindel_altopt.py -a {align_file} -n {tree_file} -o {output_folder} -m {opt_option} > {log_file}"\
                .format(align_file=align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder,opt_option='N')
    print(cmd)
    subprocess.run(cmd,shell=True)

    ## evaluation
    eval_log_file = output_folder + 'mip_evaluate_log_1.txt'
    cmd = "python metrics_mip.py -f {pr_folder} -a {align_file} -n {tree_file} > {eval_log_file}".\
        format(align_file = align_file,tree_file=tree_file,eval_log_file=eval_log_file,pr_folder = output_folder)
    print(cmd)
    subprocess.run(cmd,shell=True)

    ## Move the mip indel file to backup
    shutil.move(data_folder + 'mip_ancestor_indel.fasta', 
                data_folder + 'mip_ancestor_indel1.fasta')
    shutil.move(data_folder + 'mip_objscore.csv', 
                data_folder + 'mip_objscore1.csv')
    shutil.move(data_folder + 'mip_indscore.csv', 
                data_folder + 'mip_indscore1.csv')
    shutil.move(data_folder + 'mip_out_dist_percent.csv', 
                data_folder + 'mip_out_dist_percent1.csv')
    shutil.move(data_folder + 'mip_percent_ancestors_with_3_mut.csv', 
                data_folder + 'mip_percent_ancestors_with_3_mut1.csv')

    # MIP - 2 run
    print(f"Running MIP - SECOND INDEL SOLUTION")
    log_file  = output_folder + 'mip_inference_log_2.txt'
    cmd = "python run_mipindel_altopt.py -a {align_file} -n {tree_file} -o {output_folder} -m {opt_option} > {log_file}"\
            .format(align_file=align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder,opt_option='Y')
    print(cmd)
    subprocess.run(cmd,shell=True)

    ## evaluation
    eval_log_file = output_folder + 'mip_evaluate_log_2.txt'
    cmd = "python metrics_mip.py -f {pr_folder} -a {align_file} -n {tree_file} > {eval_log_file}".\
        format(align_file = align_file,tree_file=tree_file,eval_log_file=eval_log_file,pr_folder = output_folder)
    print(cmd)    
    subprocess.run(cmd,shell=True)

    shutil.move(data_folder + 'mip_ancestor_indel.fasta', 
                data_folder + 'mip_ancestor_indel2.fasta')
    shutil.move(data_folder + 'mip_objscore.csv', 
                data_folder + 'mip_objscore2.csv')
    shutil.move(data_folder + 'mip_indscore.csv', 
                data_folder + 'mip_indscore2.csv')
    shutil.move(data_folder + 'mip_out_dist_percent.csv', 
                data_folder + 'mip_out_dist_percent2.csv')
    shutil.move(data_folder + 'mip_percent_ancestors_with_3_mut.csv', 
                data_folder + 'mip_percent_ancestors_with_3_mut2.csv')


if __name__ == "__main__":
    eval_alt_indel()
    print("Completed")

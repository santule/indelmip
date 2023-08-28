import os
import subprocess
import sys
import getopt

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python evaluation_scripts.py -b y -m y -p y -s y')
    exit()

# main function
def main():
    run_mip = None
    run_bep = None
    run_sicp  = None
    run_psp = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "b:m:p:s:")

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

    if run_bep is None or run_mip is None or run_psp is None or run_sicp is None:
        help()

    script_folder = '/media/WorkingSpace/Share/mipindel/scripts/'
    # Run evaluation scripts

    for pr in ['CYP2U_165']:
        #,'MBL_243','CYP2U_359','GDH-#GOx_399','DHAD_585','CYP2U_595','KARI_716','KARI_1176','DHAD_1612','DHAD_1658','ALS_1990','MBL_624','CYP_3000']:

      print(f"Processing protein family {pr}")
      align_file = pr + '.aln'
      tree_file  = pr + '.nwk'
      output_folder = '.'

      print("Change directory")
      os.chdir('./'+ pr)
      print(os.getcwd())

      if run_bep == 'y':
          if pr not in ['CYP_3000','MBL_624']:
              # BEP
              ## Inference
              print(f"Processing BEP Inference")
              log_file   = 'bep_inference_log.txt'
              cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method BEP --onlyindel --save-all --threads 1 --orphans -pre 'bep' -o {output_folder}> {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
              print(cmd)
              subprocess.run(cmd,shell=True)

              log_file   = 'bep_inference_full_log.txt'
              cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method BEP --orphans -pre 'bep' -o {output_folder} > {log_file}".format(\
                        align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
              print(cmd)
              subprocess.run(cmd,shell=True)

              ## evaluation
              grasp_tree = 'bep_ancestors.nwk'
              eval_log_file = 'bep_evaluate_log.txt'
              grasp_ancestor_file = 'bep_ancestors.fa'
              cmd = "python {script_folder}grasp_evaluate.py -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"bep\" > {eval_log_file}".format(script_folder = script_folder, grasp_ancestor_file=grasp_ancestor_file,align_file = align_file,grasp_tree=grasp_tree,eval_log_file=eval_log_file)
              print(cmd)
              subprocess.run(cmd,shell=True)

      if run_psp == 'y':
          # PSP
          ## Inference
          print(f"Processing PSP Inference")
          log_file   = 'psp_inference_log.txt'
          cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method PSP --onlyindel --save-all -pre 'psp' --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
          print(cmd)
          subprocess.run(cmd,shell=True)

          log_file   = 'psp_inference_full_log.txt'
          cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method PSP -pre 'psp' --orphans -o {output_folder} > {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
          print(cmd)
          subprocess.run(cmd,shell=True)

          ## evaluation
          grasp_tree = 'psp_ancestors.nwk'
          eval_log_file = 'psp_evaluate_log.txt'
          grasp_ancestor_file = 'psp_ancestors.fa'
          cmd = "python {script_folder}grasp_evaluate.py -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"psp\" > {eval_log_file}".format(script_folder = script_folder, grasp_ancestor_file=grasp_ancestor_file,align_file = align_file,grasp_tree=grasp_tree,eval_log_file=eval_log_file)
          print(cmd)
          subprocess.run(cmd,shell=True)


      if run_sicp == 'y':
          # SICP
          ## Inference
          print(f"Processing SICP Inference")
          log_file   = 'sicp_inference_log.txt'
          cmd = "grasp -a {align_file} -n {tree_file} --time --verbose --indel-method SICP --onlyindel --save-all -pre 'sicp'  --threads 1 --orphans -o {output_folder} > {log_file}".format(\
                  align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
          print(cmd)
          subprocess.run(cmd,shell=True)

          log_file   = 'sicp_inference_full_log.txt'
          cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method SICP -pre 'sicp' --orphans -o {output_folder} > {log_file}".format(\
                    align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
          print(cmd)
          subprocess.run(cmd,shell=True)

          ## evaluation
          grasp_tree = 'sicp_ancestors.nwk'
          eval_log_file = 'sicp_evaluate_log.txt'
          grasp_ancestor_file = 'sicp_ancestors.fa'
          cmd = "python {script_folder}grasp_evaluate.py -a {grasp_ancestor_file} -e {align_file} -n {grasp_tree} -m \"sicp\" > {eval_log_file}".format(script_folder = script_folder, grasp_ancestor_file=grasp_ancestor_file,align_file = align_file,grasp_tree=grasp_tree,eval_log_file=eval_log_file)
          print(cmd)
          subprocess.run(cmd,shell=True)


      if run_mip == 'y':
          if pr not in ['CYP_3000']:
              # MIP
              ## Inference
              print(f"Running MIP")
              tree_file = 'psp_ancestors.nwk'
              log_file = 'mip_inference_log.txt'
              cmd = "python /media/WorkingSpace/Share/mipindel/scripts/main_mip_run.py -a {align_file} -n {tree_file} -o '.' > {log_file}".format(align_file=align_file,tree_file=tree_file,log_file=log_file)
              print(cmd)
              subprocess.run(cmd,shell=True)

              ## infer ancestors for MIP
              # log_file   = 'log_grasp_sicp_full.txt'
              # cmd = "grasp -a {align_file} -n {tree_file} --time --verbose -s LG --save-as FASTA --indel-method SICP -pre 'sicp' --orphans -o {output_folder} > {log_file}".format(\
              #           align_file = align_file,tree_file=tree_file,log_file=log_file,output_folder=output_folder)
              # print(cmd)
              # subprocess.run(cmd,shell=True)

              ## evaluation
              tree_file = 'psp_ancestors.nwk'
              eval_log_file = 'mip_evaluate_log.txt'
              cmd = "python {script_folder}mip_evaluate.py -e {align_file} -n {tree_file} > {eval_log_file}".format(script_folder = script_folder, align_file = align_file,tree_file=tree_file,eval_log_file=eval_log_file)
              print(cmd)
              subprocess.run(cmd,shell=True)

      os.chdir('../')
      print(os.getcwd())


if __name__ == "__main__":
    main()

import os
import subprocess
import sys
import travis_post_processing
import travis_check

def eff_pos_100(base_folder):

    d = 0.2
    s = 'A'*3554

    print("Generating synthetic indel data using Travis for 100% effective positions")
    # 100 % effective positions
    for num_extant in [200,400,600,800,1000,2000]:
        for seq_len in [100,400,600,1000]:
            for eff_pos in [100]:
                
                print(f"Travis making indel patterns with {num_extant} extants of {seq_len} each and {eff_pos} effective positions")
                isExist = os.path.exists(str(base_folder) + 't' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                if not isExist:
                    #print("Make directory")
                    os.mkdir(str(base_folder)  +  't' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                
                os.chdir(str(base_folder)  +  't' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                #print(os.getcwd())

                cmd = "travis {seq} -out all_sequences.aln -nwk input_tree.nwk -model LG  -gap -format FASTA  -seed 100 -extants {num_extants}  -dist {dist}"\
                    .format(seq = s,num_extants=num_extant,dist=d)

                #print(cmd)
                subprocess.run(cmd,shell=True)
                print("Finished")

                os.chdir('../')
                print(os.getcwd())


def eff_pos_60(base_folder):
    print("Generating synthetic indel data using Travis for 60% effective positions")

    dist_ref = {200:0.013,400:0.009,600:0.007,800:0.004,1000:0.005,2000:0.002}
    s = 'A'*3554

    # 60 % effective positions
    for num_extant in [200,400,600,800,1000,2000]:
        d = dist_ref[num_extant]
        for seq_len in [100,400,600,1000]:
            for eff_pos in [60]:
                
                print(f"Travis making indel patterns with {num_extant} extants of {seq_len} each and {eff_pos} effective positions")
                isExist = os.path.exists('./t' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                if not isExist:
                    #print("Make directory")
                    os.mkdir('./t' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                    os.chdir('./t' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                    print(os.getcwd())

                cmd = "travis {seq} -out all_sequences.aln -nwk input_tree.nwk -model LG  -gap -format FASTA  -seed 100 -extants {num_extants}  -dist {dist}"\
                    .format(seq=s,num_extants=num_extant,dist=d)

                #print(cmd)
                subprocess.run(cmd,shell=True)
                print("Finished")

                os.chdir('../')
                print(os.getcwd())


def eff_pos_20(base_folder):
    print("Generating synthetic indel data using Travis for 20% effective positions")
    s = 'A'*3554

    dist_ref = {200:0.003,400:0.0020,600:0.0015,800:0.0011,1000:0.0009,2000:0.0006}
    # 60 % effective positions
    for num_extant in [200,400,800,600,1000,2000]:
        d = dist_ref[num_extant]
        for seq_len in [100,400,600,1000]:
            for eff_pos in [20]:
                
                print(f"Travis making indel patterns with {num_extant} extants of {seq_len} each and {eff_pos} effective positions")
                isExist = os.path.exists('./t' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                if not isExist:
                    #print("Make directory")
                    os.mkdir('./t' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                    os.chdir('./t' + str(num_extant) + 'l' + str(seq_len) + 'e' + str(eff_pos))
                    print(os.getcwd())

                cmd = "travis {seq} -out all_sequences.aln -nwk input_tree.nwk -model LG  -gap -format FASTA  -seed 100 -extants {num_extants}  -dist {dist}"\
                    .format(seq=s,num_extants=num_extant,dist=d)

                #print(cmd)
                subprocess.run(cmd,shell=True)
                print("Finished")

                os.chdir('../')
                print(os.getcwd())


if __name__ == "__main__":
  base_folder = sys.argv[1] #'/Users/sanjanatule/Documents/uq/Projects/Indels/indelmip/data/travis/'
  eff_pos_100(base_folder)
  eff_pos_60(base_folder)
  eff_pos_20(base_folder)
  print("Finished running Travis")
  
  # post processing of travis
  print("Starting post processing of Travis data")
  travis_post_processing.make_extants(base_folder)

  # check if the desired dataset satisfy the requirements.
  travis_check.check_stats(base_folder,detail_print=True)

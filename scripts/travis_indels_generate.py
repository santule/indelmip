import os
import subprocess
import sys

def generate_synthetic_data(base_folder):
    s_lkp = {300:'A' * 100,500:'A' * 80,700:'A' * 60,1000:'A' * 30,2000:'A' * 20}
    for num_extant in [2000]: # total extants to create 300,500,700,1000,3
        for d in [0.02,0.5,0.8,1]: # distance parameter
            for sh in [0.5,1,2]: # shape parameter
                
                if num_extant >= 1000 and d == 0.02:
                    s = 'A' * 80
                else:
                    s = s_lkp[num_extant]
                print(f"Travis making indel patterns with {num_extant} extants of {d} distance and {sh} shape")
                isExist = os.path.exists(str(base_folder) + 't' + str(num_extant) + 'd' + str(d) + 's' + str(sh))
                if not isExist:
                    #print("Make directory")
                    os.mkdir(str(base_folder) + 't' + str(num_extant) + 'd' + str(d) + 's' + str(sh))
                
                os.chdir(str(base_folder) + 't' + str(num_extant) + 'd' + str(d) + 's' + str(sh))
                #print(os.getcwd())

                cmd = "travis {seq} -out all_sequences.aln -nwk input_tree.nwk -model LG  -gap -format FASTA  -seed 100 -extants {num_extants}  -dist {dist}\
                        -shape {sh}".format(seq = s,num_extants=num_extant,dist=d,sh=sh)

                #print(cmd)
                subprocess.run(cmd,shell=True)
                print("Finished")

                os.chdir('../')
                print(os.getcwd())

if __name__ == "__main__":
  base_folder = sys.argv[1] #'/Users/sanjanatule/Documents/uq/Projects/Indels/indelmip/data/travis/'
  generate_synthetic_data(base_folder)
  print("Finished running Travis")

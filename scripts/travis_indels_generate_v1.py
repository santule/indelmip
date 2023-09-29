import os
import subprocess
import sys
import travis_post_processing
import travis_check

def generate_synthetic_data(base_folder):
    s = 'A'*300
    for num_extant in [3000]: # 200,500,1000,3000 - total extants to create
        for d in [0.01,0.02,0.2,0.5]: # distance parameter
            for sh in [0.05,0.5,1,2]: # shape parameter
                for sc in [0.1,0.2,0.5,1,2,4]: # scale parameter
                
                    print(f"Travis making indel patterns with {num_extant} extants of {d} distance and {sh} shape and {sc} scale")
                    isExist = os.path.exists(str(base_folder) + 't' + str(num_extant) + 'd' + str(d) + 's' + str(sh) + 'c' + str(sc))
                    if not isExist:
                        #print("Make directory")
                        os.mkdir(str(base_folder) + 't' + str(num_extant) + 'd' + str(d) + 's' + str(sh)  + 'c' + str(sc))
                    
                    os.chdir(str(base_folder) + 't' + str(num_extant) + 'd' + str(d) + 's' + str(sh) + 'c' + str(sc))
                    #print(os.getcwd())

                    cmd = "travis {seq} -out all_sequences.aln -nwk input_tree.nwk -model LG  -gap -format FASTA  -seed 100 -extants {num_extants}  -dist {dist}\
                          -shape {sh} -scale {sc}".format(seq = s,num_extants=num_extant,dist=d,sh=sh,sc=sc)

                    #print(cmd)
                    subprocess.run(cmd,shell=True)
                    print("Finished")

                    os.chdir('../')
                    print(os.getcwd())



if __name__ == "__main__":
  base_folder = sys.argv[1] #'/Users/sanjanatule/Documents/uq/Projects/Indels/indelmip/data/travis/'
  generate_synthetic_data(base_folder)
  print("Finished running Travis")

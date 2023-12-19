''' main python file to run mip inference and get solution, visuals and other stats '''
import sys
import getopt
import grasp_output_process,objective_scoring,indel_scoring,new_old_patterns,check_distribution,ancestor_3_mutation_away,ancestor_2_mutation_away
from utils import write_to_file
import site_3_mut_away
from os import walk
import os

# main function
def main():
    
    prefix_method = 'tr'
    nwk_file_path = 'input_tree.nwk'
    extant_file   = 'extants.aln'
    combined_file = 'all_sequences.aln'
    combined_file_binary = 'tr_grasp_all_indel.fasta'

    # run the script from travis folder...
    for (sub_folder, _, _) in walk('.'):
        if sub_folder != '.':
            sub_folder = sub_folder

            print(f"Processing synthetic protein family {sub_folder}")

            print("Change directory")
            os.chdir(sub_folder)
            print(os.getcwd())

            # convert sequences to binary
            print("1 - COMBINE THE FASTA FILES AND CONVERT TO BINARY")
            dict_parameter = {'output_folder_location':'.','fasta_file_name':combined_file,\
                                    'grasp_indel_method':prefix_method}
            grasp_output_process.convert_grasp_fasta_to_binary(**dict_parameter)

            # Objective Scoring
            print("2 - RUNNING OBJECTIVE SCORING")
            objscore = objective_scoring.main(nwk_file_path,combined_file_binary)

            # Indel Scoring
            print("3 - RUNNING INDEL SCORING")
            indscore = indel_scoring.main(nwk_file_path,combined_file_binary)

            # New Old pattern count
            print("4 - RUNNING INDEL OLD/NEW PATTERNS")
            extant_pattern,new_pattern = new_old_patterns.main(nwk_file_path,combined_file_binary)

            # Out of distribution
            print("5 - CHECK OUT OF DISTRIBUTION PATTERNS")
            out_dist_percent = check_distribution.main(nwk_file_path,combined_file_binary,extant_file,True)

            # 3 mutation away ancestors
            print("6 - ANCESTORS 3 MUTATIONS AWAY")
            total_ancestors_i3_mut,percent_ancestors_with_3_mut = ancestor_3_mutation_away.main(nwk_file_path,combined_file_binary)

            # 2 mutation away ancestors
            print("7 - ANCESTORS 2 MUTATIONS AWAY")
            total_ancestors_i2_mut,percent_ancestors_with_2_mut = ancestor_2_mutation_away.main(nwk_file_path,combined_file_binary)

            # total sites with 3 mutation away
            print("8 - SITES 3 MUTATIONS AWAY")
            percent_total_sites_3_mut,total_sites_3_mut = site_3_mut_away.main(nwk_file_path,combined_file_binary)

            # write all numbers to the files
            print("8 - WRITING METRICS TO THE FILE")
            write_to_file(prefix_method + '_objscore.csv',objscore,prefix_method)
            write_to_file(prefix_method + '_indscore.csv',indscore,prefix_method)
            write_to_file(prefix_method + '_new_pattern.csv',new_pattern,prefix_method)
            write_to_file(prefix_method + '_out_dist_percent.csv',out_dist_percent,prefix_method)
            write_to_file(prefix_method + '_percent_ancestors_with_3_mut.csv',percent_ancestors_with_3_mut,prefix_method)
            write_to_file(prefix_method + '_percent_ancestors_with_2_mut.csv',percent_ancestors_with_2_mut,prefix_method)
            write_to_file(prefix_method + '_sites_with_3_mut.csv',total_sites_3_mut,prefix_method)

            # change directory
            os.chdir('../')
            print(os.getcwd())


if __name__ == "__main__":
  main()

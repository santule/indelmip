''' main python file to run mip inference and get solution, visuals and other stats '''
import sys
import getopt
import grasp_output_process,objective_scoring,indel_scoring,new_old_patterns,check_distribution,ancestor_3_mutation_away,ancestor_2_mutation_away
from utils import write_to_file
import site_3_mut_away

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python grasp_evaluate.py -a grasp output ancestor file | -e grasp input extant file | -n newick tree file | -m grasp indel method')
    exit()

# main function
def main():
    ancestor_file = None
    extant_file = None
    nwk_file_path = None
    grasp_indel_method  = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:e:n:m:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-a']:
            ancestor_file = arg
        elif opt in ['-e']:
            extant_file = arg
        elif opt in ['-n']:
            nwk_file_path = arg
        elif opt in ['-m']:
            grasp_indel_method = arg

    if ancestor_file is None or nwk_file_path is None or extant_file is None or grasp_indel_method is None:
        help()


    # Combine ancestor and extant sequences
    print("1 - COMBINE THE FASTA FILES AND CONVERT TO BINARY")
    grasp_output_process.main(ancestor_file,extant_file,grasp_indel_method)

    # Objective Scoring
    print("2 - RUNNING OBJECTIVE SCORING")
    objscore = objective_scoring.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')

    # Indel Scoring
    print("3 - RUNNING INDEL SCORING")
    indscore = indel_scoring.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')

    # New Old pattern count
    print("4 - RUNNING INDEL OLD/NEW PATTERNS")
    extant_pattern,new_pattern = new_old_patterns.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')

    # Out of distribution
    print("5 - CHECK OUT OF DISTRIBUTION PATTERNS")
    out_dist_percent = check_distribution.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta',extant_file,True)

    # 3 mutation away ancestors
    print("6 - ANCESTORS 3 MUTATIONS AWAY")
    total_ancestors_i3_mut,percent_ancestors_with_3_mut = ancestor_3_mutation_away.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')

    # 2 mutation away ancestors
    print("7 - ANCESTORS 3 MUTATIONS AWAY")
    total_ancestors_i2_mut,percent_ancestors_with_2_mut = ancestor_2_mutation_away.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')

    # total sites with 3 mutation away
    percent_total_sites_3_mut,total_sites_3_mut = site_3_mut_away.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')

    # write all numbers to the files
    print("8 - WRITING METRICS TO THE FILE")
    write_to_file(grasp_indel_method + '_objscore.csv',objscore,grasp_indel_method)
    write_to_file(grasp_indel_method + '_indscore.csv',indscore,grasp_indel_method)
    write_to_file(grasp_indel_method + '_new_pattern.csv',new_pattern,grasp_indel_method)
    write_to_file(grasp_indel_method + '_out_dist_percent.csv',out_dist_percent,grasp_indel_method)
    write_to_file(grasp_indel_method + '_percent_ancestors_with_3_mut.csv',percent_ancestors_with_3_mut,grasp_indel_method)
    write_to_file(grasp_indel_method + '_percent_ancestors_with_2_mut.csv',percent_ancestors_with_2_mut,grasp_indel_method)
    write_to_file(grasp_indel_method + '_sites_with_3_mut.csv',percent_total_sites_3_mut,grasp_indel_method)


if __name__ == "__main__":
  main()

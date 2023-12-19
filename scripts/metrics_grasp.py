''' main python file to run mip inference and get solution, visuals and other stats '''
import sys
import getopt
import utils_grasp,score_objective,score_indel,score_distribution,score_ancestor_3_mutation_away
from utils import write_to_file

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python grasp_evaluate.py -f data_folder | -a grasp output ancestor file | -e grasp input extant file | -n newick tree file | -m grasp indel method')
    exit()

# main function
def main():
    ancestor_file = None
    extant_file = None
    nwk_file_path = None
    grasp_indel_method  = None
    data_folder = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "f:a:e:n:m:")

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
        elif opt in ['-f']:
            data_folder = arg

    if ancestor_file is None or nwk_file_path is None or extant_file is None or grasp_indel_method is None or data_folder is None:
        help()

    # Combine ancestor and extant sequences
    print("1 - COMBINE THE FASTA FILES AND CONVERT TO BINARY")
    utils_grasp.main(ancestor_file,extant_file,grasp_indel_method,data_folder)

    # Objective Scoring
    print("2 - RUNNING OBJECTIVE SCORING")
    objscore = score_objective.score_fn(nwk_file_path,data_folder + grasp_indel_method + '_grasp_all_indel.fasta')

    # Indel Scoring
    print("3 - RUNNING INDEL SCORING")
    indscore = score_indel.score_fn(nwk_file_path,data_folder + grasp_indel_method + '_grasp_all_indel.fasta')

    # Out of distribution
    print("4 - CHECK OUT OF DISTRIBUTION PATTERNS")
    out_dist_percent = score_distribution.score_fn(nwk_file_path,data_folder + grasp_indel_method + '_grasp_all_indel.fasta',extant_file,True)

    # 3 mutation away ancestors
    print("5 - ANCESTORS 3 MUTATIONS AWAY")
    total_ancestors_i3_mut,percent_ancestors_with_3_mut = score_ancestor_3_mutation_away.score_fn(nwk_file_path,data_folder + grasp_indel_method + '_grasp_all_indel.fasta')

    # write all numbers to the files
    print("6 - WRITING METRICS TO THE FILE")
    write_to_file(data_folder + grasp_indel_method + '_objscore.csv',objscore,grasp_indel_method)
    write_to_file(data_folder + grasp_indel_method + '_indscore.csv',indscore,grasp_indel_method)
    write_to_file(data_folder + grasp_indel_method + '_out_dist_percent.csv',out_dist_percent,grasp_indel_method)
    write_to_file(data_folder + grasp_indel_method + '_percent_ancestors_with_3_mut.csv',percent_ancestors_with_3_mut,grasp_indel_method)

if __name__ == "__main__":
  main()

''' python script to evaluate mip solution '''

import sys
import getopt
import score_objective,score_indel,score_distribution,score_ancestor_3_mutation_away
from utils import write_to_file

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python metrics_mip.py -f data_folder | -a input extant file | -n newick tree file')
    exit()

# main function
def main():
    ancestor_file = None
    extant_file = None
    nwk_file_path = None
    data_folder = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "f:a:n:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-a']:
            alignment_file = arg
        elif opt in ['-n']:
            nwk_file_path = arg
        elif opt in ['-f']:
            data_folder = arg

    if nwk_file_path is None or alignment_file is None or data_folder is None:
        help()
    
    # Objective Scoring
    print("2 - RUNNING OBJECTIVE SCORING")
    objscore = score_objective.score_fn(nwk_file_path,data_folder + 'mip_ancestor_indel.fasta')

    # Indel Scoring
    print("3 - RUNNING INDEL SCORING")
    indscore = score_indel.score_fn(nwk_file_path,data_folder + 'mip_ancestor_indel.fasta')

    # Out of distribution
    print("4 - CHECK OUT OF DISTRIBUTION PATTERNS")
    out_dist_percent = score_distribution.score_fn(nwk_file_path,data_folder + 'mip_ancestor_indel.fasta',alignment_file,False)

    # 3 mutation away ancestors
    print("5 - ANCESTORS 3 MUTATIONS AWAY")
    total_ancestors_i3_mut,percent_ancestors_with_3_mut = score_ancestor_3_mutation_away.score_fn(nwk_file_path,data_folder + 'mip_ancestor_indel.fasta')

    # write all numbers to the files
    print("6 - WRITING METRICS TO THE FILE")
    write_to_file(data_folder + 'mip_objscore.csv',objscore,'mip')
    write_to_file(data_folder + 'mip_indscore.csv',indscore,'mip')
    write_to_file(data_folder + 'mip_out_dist_percent.csv',out_dist_percent,'mip')
    write_to_file(data_folder + 'mip_percent_ancestors_with_3_mut.csv',percent_ancestors_with_3_mut,'mip')

if __name__ == "__main__":
  main()

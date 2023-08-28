''' python script to evaluate mip solution '''


import sys
import getopt
import indel_scoring
import new_old_patterns
import indel_events_count
import check_distribution,ancestor_3_mutation_away
from utils import write_to_file

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python mip_evaluate.py -e input extant file | -n newick tree file')
    exit()

# main function
def main():
    ancestor_file = None
    extant_file = None
    nwk_file_path = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "e:n:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-e']:
            alignment_file = arg
        elif opt in ['-n']:
            nwk_file_path = arg

    if nwk_file_path is None or alignment_file is None:
        help()


    # Indel Scoring
    print("2 - RUNNING INDEL SCORING")
    score = indel_scoring.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # Indel event count
    print("3 - RUNNING INDEL COUNTS")
    indel_events_count.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # New Old pattern count
    print("4 - RUNNING INDEL OLD/NEW PATTERNS")
    extant_pattern,new_pattern = new_old_patterns.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # Out of distribution
    print("5 - CHECK OUT OF DISTRIBUTION PATTERNS")
    out_dist_percent = check_distribution.main(nwk_file_path,'mip_ancestor_indel.fasta',alignment_file,False)

    # 3 mutation away ancestors
    print("6 - ANCESTORS 3 MUTATIONS AWAY")
    total_ancestors_i3_mut,percent_ancestors_with_3_mut = ancestor_3_mutation_away.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # write all numbers to the files
    print("7 - WRITING METRICS TO THE FILE")
    write_to_file('mip_score.csv',score,'mip')
    write_to_file('mip_new_pattern.csv',new_pattern,'mip')
    write_to_file('mip_out_dist_percent.csv',out_dist_percent,'mip')
    write_to_file('mip_percent_ancestors_with_3_mut.csv',percent_ancestors_with_3_mut,'mip')



if __name__ == "__main__":
  main()

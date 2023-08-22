''' main python file to run mip inference and get solution, and other stats '''

import sys
import getopt
import indel_scoring
import new_old_patterns
import indel_events_count
import mipindel
import check_distribution,ancestor_3_mutation_away

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python main_mip_run.py -a alignment file | -n newick tree file | -o output file location')
    exit()

# main function
def main():
    alignment_file = None
    nwk_file_path = None
    output_file_location  = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:n:o:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-a']:
            alignment_file = arg
        elif opt in ['-n']:
            nwk_file_path = arg
        elif opt in ['-o']:
            output_file_location = arg

    if alignment_file is None or nwk_file_path is None or output_file_location is None:
        help()

    # Run MIP
    print("1 - RUNNING MIP")
    mipindel.main(alignment_file,nwk_file_path,output_file_location)

    # Indel Scoring
    print("2 - RUNNING INDEL SCORING")
    indel_scoring.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # Indel event count
    print("3 - RUNNING INDEL COUNTS")
    indel_events_count.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # New Old pattern count
    print("4 - RUNNING INDEL OLD/NEW PATTERNS")
    new_old_patterns.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # Out of distribution
    print("5 - CHECK OUT OF DISTRIBUTION PATTERNS")
    check_distribution.main(nwk_file_path,'mip_ancestor_indel.fasta',alignment_file,False)

    # 3 mutation away ancestors
    print("6 - ANCESTORS 3 MUTATIONS AWAY")
    ancestor_3_mutation_away.main(nwk_file_path,'mip_ancestor_indel.fasta')

if __name__ == "__main__":
  main()

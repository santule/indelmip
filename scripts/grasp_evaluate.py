''' main python file to run mip inference and get solution, visuals and other stats '''
import sys
import getopt
import indel_scoring
import new_old_patterns
import indel_events_count
import grasp_output_process

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

    # Indel Scoring
    print("2 - RUNNING INDEL SCORING")
    indel_scoring.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')

    # Indel event count
    print("3 - RUNNING INDEL COUNTS")
    indel_events_count.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')

    # New Old pattern count
    print("4 - RUNNING INDEL OLD/NEW PATTERNS")
    new_old_patterns.main(nwk_file_path,grasp_indel_method + '_grasp_all_indel.fasta')


if __name__ == "__main__":
  main()

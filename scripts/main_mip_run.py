''' main python file to run mip inference and get solution, visuals and other stats '''
import sys
import getopt
import indel_scoring
import new_old_patterns
import indel_events_count
import visualise_sol
import mipindel

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
    print("RUNNING MIP")
    mipindel.main(alignment_file,nwk_file_path,output_file_location)

    # Indel Scoring
    print("RUNNING INDEL SCORING")
    indel_scoring.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # Indel event count
    print("RUNNING INDEL COUNTS")
    indel_events_count.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # New Old pattern count
    print("RUNNING INDEL OLD/NEW PATTERNS")
    new_old_patterns.main(nwk_file_path,'mip_ancestor_indel.fasta')

    # Visualise solution
    print("CREATING VISUALISATION")
    visualise_sol.main(nwk_file_path,'mip_ancestor_indel.fasta',"mip")

if __name__ == "__main__":
  main()

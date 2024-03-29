''' main python file to run mip inference and get solution, and other stats '''

import sys
import getopt
import mipindel_altopt

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python run_mipindel_altopt.py -a alignment file | -n newick tree file | -o output file location | -m alternate optimal sol')
    exit()

# main function
def main():
    alignment_file = None
    nwk_file_path = None
    output_file_location  = None

    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:n:o:m:")

    if len(opts) == 0:
        help()

    for opt, arg in opts:
        if opt in ['-a']:
            alignment_file = arg
        elif opt in ['-n']:
            nwk_file_path = arg
        elif opt in ['-o']:
            output_file_location = arg
        elif opt in ['-m']:
             optimal_constraint  = arg

    if alignment_file is None or nwk_file_path is None or output_file_location is None:
        help()

    # Run MIP
    print("1 - RUNNING MIP")
    mipindel_altopt.main(alignment_file,nwk_file_path,output_file_location,optimal_constraint)

if __name__ == "__main__":
  main()

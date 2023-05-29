from pysam import FastaFile,FastxFile
import sys


def sep_extants_ancestors(input_file):

    exfile = input_file.replace('input_all','input_extants')
    anfile = input_file.replace('input_all','input_ancestors')

    with FastxFile(input_file) as fh, open(exfile, mode='w') as fout, open(anfile, mode='w') as fout1:
        for entry in fh:
            if 'A' in entry.name:
              # create a new .aln file for extants
              fout.write(str(entry) + '\n')
            else:
              # create a new .aln file for ancestors
              fout1.write(str(entry) + '\n')

    fout.close()
    fout1.close()


def main():
    # input file and output file
    input_file = sys.argv[1]
    sep_extants_ancestors(input_file)


if __name__ == "__main__":
    main()

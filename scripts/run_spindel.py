''' srinir and practer dynamic programming
    created by : chongting '''
import sequence #from binfpy
import utils
from ete3 import Tree
from collections import defaultdict
from itertools import product
import pickle
from collections import OrderedDict
import sys
import getopt

# help function
def help():
    print("Incorrect or Incomplete command line arguments")
    print('python sr_pt_indel.py -e input extant file | -n newick tree file')
    exit()

def get_aln(seqs):
    seqs = sequence.readFastaFile(seqs)
    aln = {}
    for index in range(len(seqs)):
        seq = ''.join('0' if char == '-' else '1' for char in seqs[index].sequence)
        aln.update({seqs[index].name: seq})
    alignment = list(aln.values())
    return aln, alignment


def get_tree_data(tree):
    """create neighbor dict """
    tree_file = open(tree, "r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)

    # create neighbourhood object
    ancestor_list = []
    v_list = []
    order = []
    V_father_dict = defaultdict(list)
    for n in tree.traverse():
        order.append(n.name)
        if n.is_leaf() == False:
            ancestor_list.append(n.name)
            for c in n.children:
                V_father_dict[c.name] += [n.name]
        if n.is_root() == False:
            v_list.append(n.name)

    return V_father_dict, ancestor_list, v_list, order


def remove_identical_adjacent_columns(alignment, aln):
    removed_columns = []
    alignment = [list(seq) for seq in alignment]
    current_index = 0
    original_index = 0

    while current_index < len(alignment[0]) - 1:
        current_column = [seq[current_index] for seq in alignment]
        next_column = [seq[current_index + 1] for seq in alignment]

        # check and remove the identical columns
        if current_column == next_column:
            removed_columns.append(original_index + 1)
            for seq in alignment:
                seq.pop(current_index + 1)
        else:
            current_index += 1
        original_index += 1

    result_alignment = {list(aln.keys())[index]: ''.join(alignment[index]) for index in range(len(alignment))}
    return removed_columns, result_alignment


def save_data(seq):
    with open("sp_all_indel.fasta", 'w') as file:
        for k, v in seq.items():
            file.write(f">{k}\n{v}\n")


def IndelHistory(a, T):
    """parameter: a, a path of multiple alignment file
    T, a p-ath of tree file
    return a dictionary of indel history"""
    # input alignment and trees
    aln, alignment = get_aln(a)
    tree, ancestor_list, v_list, order = get_tree_data(T)


    # remove identical adjacent columns
    removed_columns, new_alignment = remove_identical_adjacent_columns(alignment, aln)
    print(f"total removed_columns {len(removed_columns)}")
    # print(f"new_alignment {new_alignment}")

    # set all possible combination of each posisition in each v
    all_combinations = [''.join(map(str, combination)) for combination in product([0, 1], repeat=len(ancestor_list))]
    
    opt_list = []
    opt_seq_list = []

    # for every slice valid for colummn 0, get the opt
    opt0_list = []
    for time in all_combinations:
        #print("Combi trying")
        sum_sign = 0
        for v in v_list:
            father = tree[v][0]
            if v in new_alignment:
                v_value = new_alignment[v][0]
            else:
                v_value = time[ancestor_list.index(v)]
            father_value = time[ancestor_list.index(father)]
            sign = abs(int(v_value)-int(father_value))
            sum_sign += sign
        opt0_list.append(sum_sign)
    opt0 = min(opt0_list)
    opt_list.append(opt0)
    opt0_seq = all_combinations[opt0_list.index(opt0)]
    opt_seq_list.append(opt0_seq)
    print(f"Completed position 1 indel inference")

    # from 1 to n, get the opt
    for i in range(1, len(list(new_alignment.values())[0])):
        #print(f"Evaluating position {i}")
        opti_list = []
        for time in all_combinations:
            # print(f"time {time}")
            dist = 0
            for v in v_list:
                #print(f"v {v}")
                father = tree[v][0]
                if v in new_alignment:
                    v_value = new_alignment[v][i]
                    brother_v = abs(int(new_alignment[v][i-1])-int(opt_seq_list[i-1][ancestor_list.index(father)]))
                else:
                    v_value = time[ancestor_list.index(v)]
                    brother_v = abs(int(opt_seq_list[i-1][ancestor_list.index(v)]) -
                                    int(opt_seq_list[i-1][ancestor_list.index(father)]))
                father_value = time[ancestor_list.index(father)]
                sign = abs(int(v_value) - int(father_value))
                #print(f"sign {sign}")
                different = abs(sign-brother_v)
                #print(f"different {different}")
                dist += (sign*different)
                #print(f"dist {dist}")
            sign_i = dist+opt_list[i-1]
            opti_list.append(sign_i)
        '''if i == 1:
            print(opti_list)
            print(opti_list.count(1))
            indices = [index for index, value in enumerate(opti_list) if value == 1]'''
        #print(f"opti_list {opti_list}")
        opt_i = min(opti_list)
        opt_list.append(opt_i)
        opt_seq = all_combinations[opti_list.index(opt_i)]
        opt_seq_list.append(opt_seq)

    # reset the opt-seq to full result
    result_history = [''.join(nums) for nums in zip(*opt_seq_list)]
    for position in removed_columns:
        for i in range(len(result_history)):
            result_history[i] = (result_history[i][:position] + result_history[i][position - 1]
                                 + result_history[i][position:])

    result_history = [''.join(nums) for nums in result_history]
    result_history.extend(alignment)
    ancestor_list.extend(list(aln.keys()))
    all_sequence = {k: v for k, v in zip(ancestor_list, result_history)}
    all_sequence = OrderedDict((key, all_sequence[key]) for key in order)
    return save_data(all_sequence)

# main function
def main():
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

    IndelHistory(alignment_file, nwk_file_path)

if __name__ == "__main__":
  main()

from pysam import FastaFile,FastxFile
from ete3 import Tree
import csv

folder = '/media/WorkingSpace/Share/mipindel/data/real_mip_altopt/RHYS_1263/'
nwk_file_path = folder + '/psp_ancestors.nwk'

# indel solutions
#sol1 = '/media/WorkingSpace/Share/mipindel/data/real_mip_altopt/RHYS_1263/mip_ancestor_indel1.fasta'
#sol2 = '/media/WorkingSpace/Share/mipindel/data/real_mip_altopt/RHYS_1263/mip_ancestor_indel2.fasta'
#annotation_file = '/media/WorkingSpace/Share/mipindel/data/real_mip_altopt/RHYS_1263/RHYS_1263_annot_indels.csv'


# substitution solutions
sol1 = '/media/WorkingSpace/Share/mipindel/data/real_mip_altopt/RHYS_1263/mip1_ancestors.fa'
sol2 = '/media/WorkingSpace/Share/mipindel/data/real_mip_altopt/RHYS_1263/mip2_ancestors.fa'
annotation_file = '/media/WorkingSpace/Share/mipindel/data/real_mip_altopt/RHYS_1263/RHYS_1263_annot_ancestors.csv'

tree_file = open(nwk_file_path,"r")
my_tree = tree_file.read() + ";"
tree = Tree(my_tree, format=1)
sol1_out = FastaFile(sol1)
sol2_out = FastaFile(sol2)

total_diff = 0
ancestors_diff  = []
ancestors_annot = []

for n in tree.traverse(): # level order
    if not n.is_leaf():
        m_parent_sequence_1 = sol1_out.fetch(n.name)
        m_parent_sequence_2 = sol2_out.fetch(n.name)
        
        if m_parent_sequence_1 != m_parent_sequence_2:
            total_diff += 1
            ancestors_diff.append(n.name)
            ancestors_annot.append([n.name,'*'])
            
print(f"Total ancestors different {total_diff}")
print(f"Total ancestors same {1262 - total_diff}")
print(f"Ancestor names are {ancestors_diff}")

# save in csv file
with open(annotation_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(ancestors_annot)
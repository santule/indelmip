# Optimal Phylogenetic Reconstruction of Insertion and Deletion Events
A Mixed-Integer program approach to infer indel history for protein families in phylogenetic tree. We also introduce the use of partial order graph to capture multi site dependency in molecular sequences.
The resulting indel patterns can be used for Ancestral sequence reconstuction (ASR). 

<img src="https://github.com/santule/indelmip/assets/20509836/27d8b32e-e88b-43cb-a71b-ddd09a87efd8" width="400" height="400"/> 

Inputs to the method ::
* MSA (fasta format)
* Phylogenetic tree (newick format)

### Prerequisites
* GRASP  https://bodenlab.github.io/GRASP-suite/project/graspcmd/
* Gurobi https://www.gurobi.com/

### Running the method

To run the Mixed-Integer Program for indel inference for a given extant alignment file and internal branchpoint annotated phylogenetic tree

```
python run_mipindel.py -a <alignment file> -n <phylogenetic tree> -o <folder_location> -p <alpha>

arguments:
  -a  fasta format alignment file
  -n  phylogenetic tree in newick format
  -o  folder location where the output files will be stored
  -p  alpha hyperparameter value. e.g. 2
```

To evaluate MIP solution 
```
python metrics_mip.py -f <folder_location> -e  <alignment file> -n <phylogenetic tree>

arguments:
  -f  folder location where all the files for protein family are stored
  -e  ouput of the mipmodel.py fasta file consisting of ancestors and extants
  -n  phylogenetic tree in newick format
```
To run indel inference on all methods
```
python run_indel_real_datasets.py -f <folder_location> -b <y/n> -m <y/n> -p <y/n> -s <y/n>

arguments:
  -f  folder location where all the files for protein family are stored
  -b  run BEP indel inference method
  -m  run MIP indel inference method
  -p  run PSP indel inference method
  -s  run SICP indel inference method
```

To generate synthetic indels using Travis
Travis is part of GRASP suite [[1]](#1).
```
python generate_syn_indels.py <folder_location>
```

### Example solution visualisation from MIP and other methods for RNaseZ_624 family

<img src="https://github.com/santule/indelmip/assets/20509836/9a3a5840-66bf-4882-bc55-f99863e8bc31" width="500" height="500"/> 

## References
<a id="1">[1]</a> 
Foley et. al. (2022). 
Engineering indel and substitution variants of diverse and ancient enzymes using Graphical Representation of Ancestral Sequence Predictions. 
[10.1371/journal.pcbi.1010633.](https://doi.org/10.1371/journal.pcbi.1010633)

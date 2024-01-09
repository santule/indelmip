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
```

To evaluate MIP solution 
```
python metrics_mip.py -f <folder_location> -e  <alignment file> -n <phylogenetic tree>
```
To run indel inference on all methods
```
python run_indel_real_datasets.py -f <folder_location> -b <y/n> -m <y/n> -p <y/n> -s <y/n>
```

To generate synthetic indels using Travis
```
python generate_syn_indels.py <folder_location>
```

### Example solution visualisation from MIP and other methods for RNaseZ_624 family

<img src="https://github.com/santule/indelmip/assets/20509836/9a3a5840-66bf-4882-bc55-f99863e8bc31" width="500" height="500"/> 

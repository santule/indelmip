# Optimal Phylogenetic Reconstruction of Insertion and Deletion Events
A Mixed-Integer program approach to infer indel history for protein families in phylogenetic tree. We also introduce the use of partial order graph to capture multi site dependency in molecular sequences.
The resulting indel patterns can be used for Ancestral sequence reconstuction (ASR). 

<img src="https://github.com/santule/indelmip/assets/20509836/27d8b32e-e88b-43cb-a71b-ddd09a87efd8" width="450" height="400"/> 

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
  -n  ancestors annotated phylogenetic tree in newick format
  -o  folder location where the output files will be stored
  -p  alpha hyperparameter value. e.g. 2

example:
python run_mipindel.py -a ../examples/data/extants.aln -n  ../examples/data/org_tree_annonated.nwk -o ../examples/data/ -p 2
```

To evaluate MIP Indel solution
```
python metrics_mip.py -f <folder_location> -a  <alignment file> -n <phylogenetic tree>

arguments:
  -f  folder location where all the files for protein family are stored
  -a  fasta format alignment file
  -n  ancestors annotated phylogenetic tree in newick format

example:
python metrics_mip.py -f ../examples/data/ -a  ../examples/data/extants.aln -n ../examples/data/org_tree_annonated.nwk
```
To run indel inference and evulation metrics on all indel methods for a protein family
```
python run_all_indel_methods.py -a <alignment file> -n <phylogenetic tree> -o <folder_location> -b <y/n> -m <y/n> -p <y/n> -s <y/n>

arguments:
  -a  fasta format alignment file
  -n  ancestors annotated phylogenetic tree in newick format
  -o  folder location where the output files will be stored
  -b  run BEP indel inference method
  -m  run MIP indel inference method
  -p  run PSP indel inference method
  -s  run SICP indel inference method

example:
python run_all_indel_methods.py -a ../examples/data/extants.aln -n  ../examples/data/org_tree_annonated.nwk -o ../examples/data/ -b y -m y -s y -p y
```

To generate synthetic indels using Travis
Travis is part of GRASP suite [[1]](#1).
```
python generate_syn_indels.py <folder_location>

example:
python generate_syn_indels.py /data/synthetic_data/
```

To generate synthetic data with Travis using different settings
```
# TRAVIS COMMAND LINES

Usage: asr.TrAVIS
	[<ancestor-seq>]
	[-nwk <tree-file> -out <output-file-or-dir>]
	{-model <JTT(default)|Dayhoff|LG|WAG|JC|Yang>}
	{-rates <a>}
	{-seed <random>}
	{-extants <5(default)>}
	{-dist <mean-extant-to-root>
	{-shape <1.1(default)>}
	{-scale <0.2(default)>}
	{-indel <param>}
	{-gap}
	{-format <FASTA(default)|CLUSTAL|DOT|TREE|RATES|DIR>}

# create toy trees using travis
travis AAAAAAAA -out all_sequences.aln -nwk input_tree.nwk -model LG  -gap -format FASTA  -seed 200 -extants 10  -dist 0.1
```
To run MIP Indel inferece to find multiple optimal solutions (example 2 optimal solutions here)
```
python eval_mip_altopt.py -a <alignment file> -n <phylogenetic tree> -o <folder_location>

arguments:
  -a  fasta format alignment file
  -n  ancestors annotated phylogenetic tree in newick format
  -o  folder location where the output files will be stored

example:
python eval_mip_altopt.py -o ../examples/data/ -n ../examples/data/org_tree_annonated.nwk -a ../examples/data/extants.aln
```

### Example solution visualisation from MIP and other methods for RNaseZ_624 family

<img src="https://github.com/santule/indelmip/assets/20509836/9a3a5840-66bf-4882-bc55-f99863e8bc31" width="500" height="500"/> 

## References
<a id="1">[1]</a> 
Foley et. al. (2022). 
Engineering indel and substitution variants of diverse and ancient enzymes using Graphical Representation of Ancestral Sequence Predictions. 
[10.1371/journal.pcbi.1010633.](https://doi.org/10.1371/journal.pcbi.1010633)

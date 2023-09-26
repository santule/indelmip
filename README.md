# Inferring indel events in protein families using Mixed-Integer Programming for cohesive evolutionary history


To run the Mixed-Integer Program for ancestral indel inference for a given extant alignment file and internal branchpoint annotated phylogenetic tree

```
python main_mip_run.py -a CYP2U_165.aln -n CYP2U_165_ancestors.nwk -o '.'
```

To evaluate MIP solution 
```
python mip_evaluate.py -e  CYP2U_165.aln -n CYP2U_165_ancestors.nwk
```

To generate synthetic indels using Travis
```
python travis_indels_generate.py ./data/travis/
```

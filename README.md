# Optimal Phylogenetic Reconstruction of Insertion and Deletion Events.
![image](https://github.com/santule/indelmip/assets/20509836/27d8b32e-e88b-43cb-a71b-ddd09a87efd8)


To run the Mixed-Integer Program for indel inference for a given extant alignment file and internal branchpoint annotated phylogenetic tree

```
python run_mipindel.py -a alignment file -n phylogenetic tree -o folder_location -p 2
```

To evaluate MIP solution 
```
python metrics_mip.py -f folder_location -e  alignment file -n phylogenetic tree
```
To run indel inference on all methods
```
python run_indel_real_datasets.py -f folder_location -b y -m y -p y -s y
```

To generate synthetic indels using Travis
```
python generate_syn_indels.py /data/travis/
```

## Solutions from MIP and other methods visualised for RNaseZ_624

![combined_tree_plot](https://github.com/santule/indelmip/assets/20509836/9a3a5840-66bf-4882-bc55-f99863e8bc31)

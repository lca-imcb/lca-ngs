# lca-ngs
Scripts for specific NGS processing tasks

# Usage
Clone repository
```
git clone https://github.com/lca-imcb/lca-ngs
cd lca-ngs
```
Create and activate conda environment
```
conda env create -f env.yml
source activate lca-ngs
```
Read usage messages
```
script.py -h
```

# Script descriptions

`contam_filter.py` - Given BAM alignments of same reads to target and contamination genome, 
remove contaminant reads and perform quality filtering. Contamination is detected by 
comparing mean mapping qualities of all aligned segments for each read.

`bam_split_by_readgroup.py` - Split bam file based on read groups.

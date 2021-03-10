#!/usr/bin/env python3

import sys
import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from os import path

# Set seaborn plot style
sns.set_style("white")

def transform_func(x):
    '''Transform genotype information to a numerical value'''
    if x == './.':
        return 0
    elif x == '0/0':
        return 0
    elif x == '1/0':
        return 1
    elif x == '0/1':
        return 2
    elif x == '1/1':
        return 3

# Read command line arguments
vcf_file = sys.argv[1]
output_dir = sys.argv[2]

name = path.basename(vcf_file).split('.gz')[0].split('.vcf')[0]

# Read vcf file
lines = []
with gzip.open(vcf_file, 'rb') as f:
    for line in f:
        line_decoded = line.decode()
        if '##' in line_decoded:
            continue
        else:
            lines.append(line_decoded.strip('\n').split('\t'))
vcf = pd.DataFrame(lines[1:], columns=lines[0])

# Apply transform_func to all sample columns from vcf
vcf_transformed = vcf.iloc[:,9:].applymap(lambda x: transform_func(x))

# Compute pairwise pearson correlations between all possible sample pairs
vcf_corr = vcf_transformed.corr()

# Plot the correlation matrix from above as a heatmap and save as png
plt.figure(figsize=(15, 15), dpi=150)
hmap = sns.heatmap(vcf_corr, square=True, cmap='YlGnBu', linewidth=.2, cbar_kws={"shrink": 0.5})
hmap.figure.savefig(
   path.join(output_dir, f'{name}.corrmap.png'),
   bbox_inches='tight',
   transparent=False
)
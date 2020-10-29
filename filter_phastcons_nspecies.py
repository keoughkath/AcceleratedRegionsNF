# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Filter phastCons elements for those with a minimum number of species in their multiple alignment. 

Usage:
    filter_phastcons_nspecies.py <phastcons_path> <ss_dir> <min_n_species> <phastcons_out>

Arguments:
    phastcons          phastCons elements BED file path
    ss_dir               Directory where SS-formatted MSAs for phastCons are stored
    min_n_species      Minimum # of species required in phastCons MSA for it to pass filtering
    phastcons_out       Filtered phastCons returned in BED format

"""

from docopt import docopt
import pandas as pd
import os

"""
Developed by Kathleen Keough, Pollard Lab @ Gladstone
2020
"""


"""
FUNCTION DEFINITIONS
"""

# get number of species in SS file for each one

def get_n_species(row, ssdir):
    uid = row['name']
    with open(os.path.join(ssdir, f'{uid}.ss'), 'r') as f:
        n_species = int(f.readline().split(' = ')[-1])
    return(n_species)

"""
MAIN
"""

def main():

    # specify inputs

    phastcons_path = args['<phastcons_path>']
    ss_dir = args['<ss_dir>'] 
    min_n_species = int(args['<min_n_species>'])

    # load phastCons

    phastcons = pd.read_table(phastcons_path, header=None, 
        names=['chrom','start','end','name'])

    # get # of species per phastCons MSA

    phastcons['n_species'] = phastcons.apply(lambda row: get_n_species(row, ss_dir), axis=1)

    # filter phastcons for those that have greater than specified # of species

    phastcons_filt = phastcons.query('n_species >= @min_n_species')[['chrom','start','end','name']].copy()

    # save to file

    phastcons_filt.to_csv(args['<phastcons_out>'], 
        sep='\t', index=False, header=False)
    

if __name__ == "__main__":
    args = docopt(__doc__)
    main()


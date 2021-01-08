# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Filter phastCons elements for those with a minimum number of species in their multiple alignment. 

Usage:
    filter_phastcons_nspecies.py <phastcons_in> <phastcons_out>

Arguments:
    phastcons_in          phastCons elements BED file path
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

def score_filter(indf):
	indf['size'] = indf['end'] - indf['start']
	indf['score_norm'] = indf['score'] / indf['size']
	indf['decile'] = pd.qcut(indf['score_norm'], 10, labels=False)
	return(indf.query('decile > 3.0')[['chrom','start','end','name']])

def main():
	# specify inputs
	phastcons_in = pd.read_csv(args['<phastcons_in>'], sep='\t', 
		names=['chrom','start','end','name','score','strand'])

	# filter by score
	phastcons_out = score_filter(phastcons_in)

	# save to file
	phastcons_out.to_csv(args['<phastcons_out>'], 
        sep='\t', index=False, header=False)

if __name__ == "__main__":
    args = docopt(__doc__)
    main()


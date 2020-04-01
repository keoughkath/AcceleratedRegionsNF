# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Take phastCons elements scored for acceleration by phyloP, 
perform Benjamini-Hochberg multiple test correction of p-values, 
and filter by user-specified significance threshold.

Usage:
	multiple_test_correction_filtering.py <phastcons> <max_p> <AR_out> <phastCons_out>

Arguments:
	phastcons          phastCons elements scored by phyloP for acceleration (GFF file)
	max_p			   maximum p-value allowed after multiple test correction
	AR_out             where to save the significance-filtered phastCons elements (ARs!)
	phastCons_out	   all scored phastCons returned

"""

from docopt import docopt
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')


def get_pvalue(log10p):
	# convert from log 10 p-values to normal p-values
	return 10**-log10p


def main(args):
	print(args)
	phastcons = pd.read_csv(args['<phastcons>'], sep='\t', low_memory=False, header=1, names=['chrom',
			                                                                                 'ftype',
			                                                                                 'feature_type',
			                                                                                 'start_',
			                                                                                 'stop',
			                                                                                 '-log10p',
			                                                                                 'null1',
			                                                                                 'null2',
			                                                                                 'id'])

	# convert coordinates from 1-based (GFF) to 0-based (BED)
	phastcons['start'] = phastcons['start_'] - 1

	# convert from log 10 p-values to normal p-values (unadjusted)
	phastcons['pval'] = phastcons['-log10p'].apply(get_pvalue)

	# get BH-corrected p-values
	phastcons['p_adjust'] = stats.p_adjust(FloatVector(phastcons['pval'].tolist()), method = 'BH')

	# filter for corrected p-values below user-specified threshold (these are the "official" accelerated elements)
	max_p = float(args['<max_p>'])
	ARs = phastcons.query('p_adjust < @max_p').copy()
	ARs['id'] = ARs.index
	print(f'{len(ARs)} ARs ({round(100.0*(len(ARs) / len(phastcons)), 2)}% of phastCons) based on significance threshold of {max_p}')

	# save ARs to BED file
	ARs[['chrom','start','stop','p_adjust']].to_csv(args['<AR_out>'], sep='\t', header=False, index=False)
	phastcons.to_csv(args['<phastCons_out>'], sep='\t', index=False)


if __name__ == '__main__':
	args = docopt(__doc__)
	main(args)
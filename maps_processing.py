from confidence_calculation import *
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import cooler
import bioframe
import cooltools
import cooltools.eigdecomp
import cooltools.expected


def average_tad(file_names, k, tads_positions, area_size):
	mtx = {}
	confidence_mtx = {}
	avg_tad = {}
	mtx['rep_1'], regions, bin_size = download_matrix(file_names[0])
	mtx['rep_2'], regions, bin_size = download_matrix(file_names[1])
	for chrm in regions:
		confidence_mtx[chrm] = get_confidence(mtx['rep_1'](chrm,chrm),mtx['rep_2'](chrm,chrm),k)
	for rep in ['rep_1','rep_2']:
		tads_areas = []
		confidence_areas = []
		for chrm in regions:
			# Left border
			tads_areas += [extract_tad(mtx[rep](chrm,chrm), tads_positions.loc[tad,:], bin_size, area_size, border='Left', confidence=False) for tad in tads_positions[tads_positions['Chr']==chrm[0]].index]
			confidence_areas += [extract_tad(confidence_mtx[chrm], tads_positions.loc[tad,:], bin_size, area_size, border='Left', confidence=True) for tad in tads_positions[tads_positions['Chr']==chrm[0]].index]
			# Right border
			tads_areas += [extract_tad(mtx[rep](chrm,chrm), tads_positions.loc[tad,:], bin_size, area_size, border='Right', confidence=False) for tad in tads_positions[tads_positions['Chr']==chrm[0]].index]
			confidence_areas += [extract_tad(confidence_mtx[chrm], tads_positions.loc[tad,:], bin_size, area_size, border='Right', confidence=True) for tad in tads_positions[tads_positions['Chr']==chrm[0]].index]
		tads_areas = np.asarray([i for i in tads_areas if i.shape[0] == area_size*2+1])
		confidence_areas = np.asarray([i for i in confidence_areas if i.shape[0] == area_size*2+1])
		confidence_areas /= np.nansum(confidence_areas,axis=0)
		
		# Average TAD	
		avg_tad[rep] = np.nansum(confidence_areas * tads_areas,axis=0)
	return avg_tad

def choose_k(file_names, max_k = 5, window_size = 5):
	correlations = replicates_correlation(window_size, max_k, file_names)
	k = np.array(correlations).argmax()+1
	return k

def replicates_correlation(window_size, max_k, file_names):
	ins = {}   
	mtx = {}
	mtx['rep_1'], regions, bin_size = download_matrix(file_names[0])
	mtx['rep_2'], regions, bin_size = download_matrix(file_names[1])
	cor = []
	for k in range(1,max_k):
		for rep in ['rep_1','rep_2']:
			ins[rep] = []
			for chrm in regions:
				confidence_mtx = get_confidence(mtx['rep_1'](chrm,chrm),mtx['rep_2'](chrm,chrm),k)
				norm_mtx = confidence_mtx*mtx[rep](chrm,chrm)
				ins[rep] += [calc_ins(confidence_mtx,i,window_size,norm_mtx) for i in range(window_size,norm_mtx.shape[0])]
			ins[rep] = (ins[rep] - np.nanmean(ins[rep]))/np.nanstd(ins[rep])
			ins[rep] = np.nan_to_num(ins[rep])
		cor += [spearmanr(ins['rep_1'].flatten(),ins['rep_2'].flatten())[0]]
	return cor

def downsampling(file_names, k, multiplier, assembly_v, generate_files):
	Q_LO = 0.025 
	Q_HI = 0.975 
	N_GROUPS = 11 
	q_edges = np.linspace(Q_LO, Q_HI, N_GROUPS+1)

	c = []
	c += [cooler.Cooler(file_names[0])]
	c += [cooler.Cooler(file_names[1])]
	f = []
	f += [c[0].matrix(balance=True)] 
	f += [c[1].matrix(balance=True)] 

	bins = c[0].bins()[:]
	bin_size = c[0].binsize
	genecov = bioframe.tools.frac_gene_coverage(bins, assembly_v)
	regions = [(chrom, 0, c[0].chromsizes[chrom]) for chrom in c[0].chromnames]
	
	cis_eigs = []	
	cis_eigs += [cooltools.eigdecomp.cooler_cis_eig(
			c[0],
			genecov,
			regions=None,
			n_eigs=5,
			phasing_track_col='gene_count')]
	cis_eigs += [cooltools.eigdecomp.cooler_cis_eig(
			c[1],
			genecov,
			regions=None,
			n_eigs=5,
			phasing_track_col='gene_count')]
	
	group_E1_bounds = []
	group_E1_bounds += [cooltools.saddle.quantile(cis_eigs[0][1]['E1'], q_edges)]
	group_E1_bounds += [cooltools.saddle.quantile(cis_eigs[1][1]['E1'], q_edges)]

	digitized = []
	d = cooltools.saddle.digitize_track(group_E1_bounds[0],track=(cis_eigs[0][1], 'E1'))
	digitized += [d[0]]
	d = cooltools.saddle.digitize_track(group_E1_bounds[1],track=(cis_eigs[1][1], 'E1'))
	digitized += [d[0]]

	expected = []	
	expected += [cooltools.expected.cis_expected(c[0], regions, use_dask=True)]
	expected += [cooltools.expected.cis_expected(c[1], regions, use_dask=True)]

	getmatrix = []
	getmatrix += [cooltools.saddle.make_cis_obsexp_fetcher(c[0], (expected[0], 'balanced.avg'))]
	getmatrix += [cooltools.saddle.make_cis_obsexp_fetcher(c[1], (expected[1], 'balanced.avg'))]

	confidence = {reg[0] : get_confidence(getmatrix[0](reg,reg),getmatrix[1](reg,reg),k) for reg in regions}
	
	downsample(generate_files[0],
				f[0],
				bin_size, multiplier,
				(digitized[0], 'E1' + '.d'),
				contact_type='cis',
				confidence = confidence)
	downsample(generate_files[1],
				f[1],
				bin_size, multiplier,
				(digitized[1], 'E1' + '.d'),
				contact_type='cis',
				confidence = confidence)
	


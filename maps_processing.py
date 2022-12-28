from confidence_calculation import *
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import cooler
import bioframe
import cooltools
import cooltools.eigdecomp
import cooltools.expected


def average_tad(file_names, k_parameter, tads_positions, area_size):
    ''' Function for calculation of average picture of Hi-C signal around tads boundaries 

        Parameters
        ----------
        file_names : list of str 
            Path to two Hi-C replicates
        k_parameter : int 
            Special parameter for confidence calculation
        tads_positions : pandas.DataFrame
            Positions of all tads with columns <Chr,Start,End>
        area_size : int 
            Size of a window around tads boundaries to extract 

        Returns
        -------
        avg_tad : dict of two numpy.ndarray
            Dict with 2 arrays of size area_size x area_size with averaged values from original 2 Hi-C matrices around tads boundaries'''

    # Download Hi-C replicates
    mtx = {}
    mtx['rep_1'], regions, bin_size = download_matrix(file_names[0]) 
    mtx['rep_2'], regions, bin_size = download_matrix(file_names[1]) 

    # Calculate confidence parameter for every chromosome
    confidence_chrm = {}
    for chrm in regions: 
        confidence_chrm[chrm] = calculate_confidence(mtx['rep_1'](chrm,chrm),mtx['rep_2'](chrm,chrm),k_parameter)

    # Normalization to confidence and average tad calculation
    avg_tad = {} 
    for rep in ['rep_1','rep_2']:
        tads_areas = [] # list of areas of Hi-C values around tads
        confidence_areas = [] # list of areas of confidence values around tads
        for chrm in regions: 
            # Extract areas around tads 
            # Left boundaries
            tads_areas += [extract_tad(mtx[rep](chrm,chrm), tads_positions.loc[tad,:], bin_size, area_size, boundary='Left') for tad in tads_positions[tads_positions['Chr']==chrm[0]].index]
            confidence_areas += [extract_tad(confidence_chrm[chrm], tads_positions.loc[tad,:], bin_size, area_size, boundary='Left') for tad in tads_positions[tads_positions['Chr']==chrm[0]].index]
            # Right boundaries
            tads_areas += [extract_tad(mtx[rep](chrm,chrm), tads_positions.loc[tad,:], bin_size, area_size, boundary='Right') for tad in tads_positions[tads_positions['Chr']==chrm[0]].index]
            confidence_areas += [extract_tad(confidence_chrm[chrm], tads_positions.loc[tad,:], bin_size, area_size, boundary='Right') for tad in tads_positions[tads_positions['Chr']==chrm[0]].index]

        # Generate numpy arrays and check if all areas have the proper size
        tads_areas = np.asarray([i for i in tads_areas if i.shape[0] == area_size*2+1]) 
        confidence_areas = np.asarray([i for i in confidence_areas if i.shape[0] == area_size*2+1])
        confidence_areas /= np.nansum(confidence_areas,axis=0) # Normalize confidence for every bin
		
        # Averaging through all tads with normalization to confidence
        avg_tad[rep] = np.nansum(confidence_areas * tads_areas,axis=0)

    return avg_tad

def choose_k(file_names, max_k = 5, window_size = 5):
    ''' Calculation of parameter k for confidence formulae, it takes two replicates, itteratively use different k for normalization and finds the best k based on the correlation of normalized replicates.

        Parameters
        ----------
        file_names : list of str 
            Path to two Hi-C replicates
        max_k : int 
            Maximum value for k
        window_size : int 
            Size of the window for insulation score calculation

        Returns
        -------
        k : int 
            k parameter value '''

    correlations = replicates_correlation(window_size, max_k, file_names) # The list of replicates' correlations for every k used

    k = np.array(correlations).argmax()+1 # Which k gives highest replicates' correlation

    return k

def downsampling(file_names, k_parameter, factor, output_file_name):
    ''' Function to downscale the Hi-C map taking into acount confidence parameter 
  
        Parameters
        ----------
        file_names : list of str 
            Path to two Hi-C replicates
        k_parameter : int 
            Special parameter for confidence calculation
        factor : 
            Coefficient for downsampling, the shape of resulting map will be equal to original shape / factor
        output_file_name : str
            Name of file to write results to

        Returns
        -------
        The result of this function is generated <file_name> file with normalized contacts for downsampled Hi-C map. 
        The resulted file is a tsv pairs file with columns: chr1 start1 end1 chr2 start2 end2 value '''

    # Download Hi-C maps for replicates
    cooler_rep1 = cooler.Cooler(file_names[0])
    cooler_rep2 = cooler.Cooler(file_names[1])
    balanced_mtx1 = cooler_rep1.matrix(balance=True) 
    balanced_mtx2 = cooler_rep2.matrix(balance=True)

    bin_size = cooler_rep1.binsize # Size of bins
    regions = [(chrom, 0, cooler_rep1.chromsizes[chrom]) for chrom in cooler_rep1.chromnames] # Names and sizes of chromosomes

    # Observed over expected
    expected_rep1 = cooltools.expected.cis_expected(cooler_rep1, regions, use_dask=True) 
    expected_rep2 = cooltools.expected.cis_expected(cooler_rep2, regions, use_dask=True)

    getmatrix_rep1 = cooltools.saddle.make_cis_obsexp_fetcher(cooler_rep1, (expected_rep1, 'balanced.avg'))
    getmatrix_rep2 = cooltools.saddle.make_cis_obsexp_fetcher(cooler_rep2, (expected_rep2, 'balanced.avg'))

    # Calculate confidence
    confidence = {chrm[0] : calculate_confidence(getmatrix_rep1(chrm,chrm),getmatrix_rep2(chrm,chrm),k_parameter) for chrm in regions}

    # Downsampling
    downsample(output_file_name+'_rep1',
               balanced_mtx1,
               bin_size, factor,
               regions,
               confidence)
    downsample(output_file_name+'_rep2',
               balanced_mtx2,
               bin_size, factor,
               regions,
               confidence)
	


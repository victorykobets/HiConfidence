import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import cooler
import bioframe
import cooltools
import cooltools.eigdecomp
import cooltools.expected
import cooltools.saddle
import warnings
warnings.filterwarnings("ignore")

def download_matrix(file_name):
    ''' Function for Hi-C matrix download 

        Parameters
        ----------
        file_name : str 
            Path to Hi-C matrix file

        Returns
        -------
        getmatrix : function
            Observed over expected Hi-C map
        regions : list of tuples of (chrm, 0, chrm_size)
            List of chromosome names and sizes
        bin_size : int
            Size of bins in matrix '''

    # Download Hi-C map
    c = cooler.Cooler(file_name)

    bin_size = c.binsize # Size of bins

    regions = [(chrom, 0, c.chromsizes[chrom]) for chrom in c.chromnames] # Names and sizes of chromosomes

    expected = cooltools.expected.cis_expected(c, regions, use_dask=True)
    getmatrix = cooltools.saddle.make_cis_obsexp_fetcher(c, (expected, 'balanced.avg')) # Observed over expected map 

    return getmatrix, regions, bin_size 	

def calculate_confidence(mtx1,mtx2,k_parameter):
    ''' Main function for calculation of confidence parameter

        Parameters
        ----------
        mtx1 : numpy.ndarray
            Hi-C map of specific chromosome, first replicate
        mtx2 : numpy.ndarray
            Hi-C mapof specific chromosome, second replicate
        k_parameter : int 
            Special parameter for confidence calculation 

        Output: 
        confidence : numpy.ndarray
            Matrix of confidence parameter calculated for every bin of two orignal Hi-C matrices, has the same shape as original matrices '''

    def diag_indices(mtx, diagonal_n):
        ''' Additional function to get indexes of bins located at specific map's diagonal '''
        rows, cols = np.diag_indices_from(mtx)
        if diagonal_n < 0:
            return rows[-diagonal_n:], cols[:diagonal_n]
        elif diagonal_n > 0:
            return rows[:-diagonal_n], cols[diagonal_n:]
        else:
            return rows, cols

    mask = np.nanmean(np.array([mtx1,mtx2]),axis=0) == 0 # To mask bins with zero values on two maps
    confidence = 1/(((np.abs(mtx2-mtx1))/(np.nanmean(np.array([mtx1,mtx2]),axis=0) + mask)+1)**k_parameter) # Confidence parameter calculation
    confidence[mask] = 0 # Masking of zero values

    # To mask two main diagonals 
    for diagonal in range(-1,2):
        confidence[diag_indices(confidence, diagonal)] = 0
    return confidence

def replicates_correlation(window_size, max_k, file_names):
    ''' Correlation of insulation scores of replicates after MuReSiE normalization 

        Parameters
        ----------
        window_size : int 
            Size of the window for insulation score calculation
        max_k : int
            Maximum value for k
        file_names : list of str 
            Path to two Hi-C replicates
 
        Returns
        -------
        correlations : list of float
            List of spearman's correlations of replicates for k from 1 to max_k
        '''

    # Download Hi-C replicates
    mtx = {} 
    mtx['rep_1'], regions, bin_size = download_matrix(file_names[0]) 
    mtx['rep_2'], regions, bin_size = download_matrix(file_names[1]) 
    
    correlations = [] # An empty list for correlations

    for k in range(1,max_k):
        # For each k from 1 to max_k calculate correlation of replicates
        insulations = {} # Dictionary for insulation scores for two replicates
        for rep in ['rep_1','rep_2']:
            insulations[rep] = []
            for chrm in regions:
                # For every chromosome MuReSiE normalization of original matrix
                confidence_mtx = calculate_confidence(mtx['rep_1'](chrm,chrm),mtx['rep_2'](chrm,chrm),k)
                norm_mtx = confidence_mtx*mtx[rep](chrm,chrm)
                # Calculation of insulation score for normalized map
                insulations[rep] += [calc_point_insulation(confidence_mtx,i,window_size,norm_mtx) for i in range(window_size,norm_mtx.shape[0])]

            insulations[rep] = (insulations[rep] - np.nanmean(insulations[rep]))/np.nanstd(insulations[rep]) # Normalization of insulation score
            insulations[rep] = np.nan_to_num(insulations[rep]) # NaNs to zeroes
        correlations += [spearmanr(insulations['rep_1'].flatten(),insulations['rep_2'].flatten())[0]] # Spearman correlation of insulation scores of two replicates

    return correlations


def calc_point_insulation(confidence_mtx,pos,window_size,hic_mtx):
    ''' Calculation of average Hi-C signal in a window near main diagonal at specific position, taking into acount confidence parameter. 

        Parameters
        ----------
        confidence_mtx : numpy.ndarray
            Confidence parameter matrix
        pos : int 
            Position
        window_size : int 
            Size of the window for the insulation score calculation
        hic_mtx : numpy.ndarray
            Original Hi-C matrix 

        Returns
        -------
        a value for an average Hi-C signal : float'''

    return np.nansum(hic_mtx[pos-window_size:pos,pos+1:pos+window_size+1]/np.nansum(confidence_mtx[pos-window_size:pos,pos+1:pos+window_size+1]))

def extract_tad(mtx, tad, bin_size, area_size, boundary = 'Left'):
    ''' Function to extract a part of Hi-C map surrounding specific tad's boundary position 

        Parameters
        ----------
        mtx : numpy.ndarray
            Hi-C matrix,
        tad : pandas.Series
            Tad's position, with columns <Chr,Start,End>
        bin_size : int
            Size of bins in matrix
        area_size : int 
            Size of a window around tad's boundary to extract
        boundary : 'Left' or 'Right'
            Left or right boundary to calculate 

        Returns
        -------
        around_tad : numpy.ndarray
            An array of size area_size x area_size with values from original Hi-C matrix around tad's boundary '''
    
    start = tad.Start//bin_size 	# Tad start position
    end = tad.End//bin_size 		# Tad end position
    shape = mtx.shape[0] 		# Shape of a Hi-C map

    around_tad = np.array([[np.nan]*(area_size*2+1)]*(area_size*2+1)) # Generate an empty array

    # Extraction of left tad's boundary
    if boundary == 'Left':
        # If left tad's boundary is too close to the left end of Hi-C map        
        if start - area_size < 0:
            around_tad[0-(start-area_size):,0-(start-area_size):] = mtx[
                              0:start+area_size+1,
                              0:start+area_size+1]
        # If area size around left tad's boundary is big enough to capture right boundary as well  
        elif start + area_size + 1 > end + 1:
            around_tad[:(end + 1 - (start + area_size + 1)),:(end + 1 - (start + area_size + 1))] = mtx[
                              start-area_size:end+1,
                              start-area_size:end+1]
        # If everything is ok
        else:
            around_tad = mtx[start-area_size:start+area_size+1,
                             start-area_size:start+area_size+1]

    # Extraction of right tad's boundary
    elif boundary == 'Right':
        # If area size around right tad's boundary is big enough to capture left boundary as well
        if end - area_size < start:
            around_tad[(start)-(end-area_size):,(start)-(end-area_size):] = mtx[
                              start:end+area_size+1,
                              start:end+area_size+1]
        # If right tad's boundary is too close to the right end of Hi-C map
        elif end + area_size + 1 > shape:
            around_tad[:(shape-(end + area_size + 1)),:(shape-(end + area_size + 1))] = mtx[
                              end-area_size:end+area_size+1,
                              end-area_size:end+area_size+1]
        # If everything is ok
        else:
            around_tad = mtx[end-area_size:end+area_size+1,
                             end-area_size:end+area_size+1]
            # Transpose of right boundary to be able to calculate average picture for left and right boundaries together
            around_tad = around_tad[::-1].transpose()[::-1] 

    return around_tad


def downsample(file_name,
    getmatrix,
    bin_size,
    factor,
    regions,
    confidence
):
    ''' Function to downscale the Hi-C map taking into acount confidence parameter. This function takes the area on the original Hi-C map with shape of factor x factor and performs an averaging taking into account the confidence of bins.

        Parameters
        ----------
        file_name : str 
            Name of file to write results to
        getmatrix : function 
            Observe over expected map fetcher
        bin_size : int 
            Size of bins in matrix
        factor : int 
            Coefficient for downsampling, the shape of resulting map will be equal to original shape / factor
        regions : list of tuples of (chrm, 0, chrm_size)
            List of chromosome names and sizes
        confidence : dict of numpy.ndarray 
            Dict of matrices of confidence parameter

        Returns
        -------
        The result of this function is generated <file_name> file with normalized contacts for downsampled Hi-C map. 
        The resulted file is a tsv pair file with columns: chr1 start1 end1 chr2 start2 end2 value '''


    f = open(file_name+'.tsv','w') 			# Open file for writing

    # For each chromosome calculate normalized values
    for chrm in regions:
        matrix = getmatrix.fetch(chrm[0]) 		# Hi-C map for specific chromosome
        confidence_ = confidence[chrm[0]].copy() 	# Confidence values for specific chromosome
        shape = matrix.shape[0] 			# shape of original Hi-C map

        # For all areas of shape factor x factor ...
        for i in range(1,shape//factor+1):
                for j in range(i,shape//factor+1):
                    # ... perform averaging taking into account confidence parameter and put the value in the file
                    if np.nansum(matrix[(i-1)*factor:i*factor,(j-1)*factor:j*factor]) > 0:
                        f.write('\t'.join([chrm[0],
                                      str((i-1)*factor*bin_size),
                                      str((i)*factor*bin_size),
                                      chrm[0],
                                      str((j-1)*factor*bin_size),
                                      str((j)*factor*bin_size),
                                      str(np.nansum(matrix[(i-1)*factor:i*factor,(j-1)*factor:j*factor]*(confidence_[(i-1)*factor:i*factor,(j-1)*factor:j*factor]/np.nansum(confidence_[(i-1)*factor:i*factor,(j-1)*factor:j*factor]))))])+'\n')

    f.close()

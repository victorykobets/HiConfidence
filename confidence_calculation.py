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
    c = cooler.Cooler(file_name)
    bins = c.bins()[:]
    bin_size = c.binsize
    regions = [(chrom, 0, c.chromsizes[chrom]) for chrom in c.chromnames]
    expected = cooltools.expected.cis_expected(c, regions, use_dask=True)
    getmatrix = cooltools.saddle.make_cis_obsexp_fetcher(c, (expected, 'balanced.avg'))
    return getmatrix, regions, bin_size

def get_confidence(df0,df1,k):    
    mask = np.nanmean(np.array([df0,df1]),axis=0) == 0
    confidence = 1/(((np.abs(df1-df0))/(np.nanmean(np.array([df0,df1]),axis=0) + mask)+1)**k)
    confidence[mask] = 0
    for k in range(-1,2):
        confidence[kth_diag_indices(confidence, k=k)] = 0
    return confidence

def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols

def calc_ins(confidence_mtx,i,window_size,mtx):
    return np.nansum(mtx[i-window_size:i,i+1:i+window_size+1]/np.nansum(confidence_mtx[i-window_size:i,i+1:i+window_size+1]))

def extract_tad(mtx, tad, bin_size, area_size, border = 'Left', confidence=False):
        start = tad.Start//bin_size
        end = tad.End//bin_size
        shape = mtx.shape[0]
        pos = np.array([[np.nan]*(area_size*2+1)]*(area_size*2+1))
        if border=='Left':
                if start-area_size < 0:
                    pos[0-(start-area_size):,0-(start-area_size):] = mtx[
                              0:start+area_size+1,
                              0:start+area_size+1]
                elif start + area_size + 1 > end + 1:
                    pos[:(end + 1 - (start + area_size + 1)),:(end + 1 - (start + area_size + 1))] = mtx[
                              start-area_size:end+1,
                              start-area_size:end+1]
                else:
                    pos = mtx[
                              start-area_size:start+area_size+1,
                              start-area_size:start+area_size+1]
        elif border == 'Right':
                if end + area_size < start:
                    pos[(start)-(end-area_size):,(start)-(end-area_size):] = mtx[
                              start:end+area_size+1,
                              start:end+area_size+1]
                elif end + area_size + 1 > shape:
                    pos[:(shape-(end + area_size + 1)),:(shape-(end + area_size + 1))] = mtx[
                              end-area_size:end+area_size+1,
                              end-area_size:end+area_size+1]
                else:
                    pos = mtx[
                              end-area_size:end+area_size+1,
                              end-area_size:end+area_size+1]
                pos = pos[::-1].transpose()[::-1]
        return pos

def downsample(file_name,
    getmatrix,
    bin_size,
    multiplier,
    digitized,
    contact_type,
    confidence,
    regions=None
):
    digitized_df, name_ = digitized
    regions = [
            (chrom, df.start.min(), df.end.max())
            for chrom, df in digitized_df.groupby("chrom")
              ]
    f = open(file_name,'w')
    n = 0
    k = 0

    for reg1 in regions:
        big_k = k
        big_n = n
        matrix = getmatrix.fetch(reg1[0])
        confidence_ = confidence[reg1[0]].copy()
        shape = matrix.shape[0]
        for i in range(1,shape//multiplier+1):
                k = big_k+i-1
                for j in range(i,shape//multiplier+1):
                    f.write('\t'.join([reg1[0],
                                      str((i-1)*multiplier*bin_size),
                                      str((i)*multiplier*bin_size),
                                      reg1[0],
                                      str((j-1)*multiplier*bin_size),
                                      str((j)*multiplier*bin_size),
                                      str(np.nansum(matrix[(i-1)*multiplier:i*multiplier,(j-1)*multiplier:j*multiplier]*(confidence_[(i-1)*multiplier:i*multiplier,(j-1)*multiplier:j*multiplier]/np.nansum(confidence_[(i-1)*multiplier:i*multiplier,(j-1)*multiplier:j*multiplier]))))])+'\n')
                    k += 1
                n += 1
    f.close()

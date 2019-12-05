"""
Code implementing algorithm from paper:
"A Peak Synchronization Measure for Multiple Signals" (Biswas et al. 2013)
https://ieeexplore.ieee.org/document/6845368/

Author: Pedro Bordignon May 2018

Modified: SP, 9/2019
"""

import sys, os
import pickle


import numpy as np
import pandas as pd

from math import floor

from scipy.stats import norm
from scipy.special import comb

import multiprocessing as mp


################################################################################
### Consensus Profile with Biswas algorithm:
################################################################################

def normcdf(x, mu, var):
    """
    Computes univariate nomal distribution function at x, corresponding to mean
    mu and variance sigma^2. (Cumulative Distribution Function)
    mu = 0
    variance = 1
    """
    return norm.cdf(x, mu, var)


def norminv(x, mu, var):
    """
    Computes inverse of univariate normal distribution function with mean u and variance var.
    (Normal inverse cumulative distribution function)
    """
    return norm.ppf(x, mu, var)


def horzcat(A, B):
    """
    Horizontally concatenates two matrices A and B to form a single matrix.
    """
    # B = np.array([1, 2]) Matrix in numpy array in python.
    return np.hstack((A, B))


def find_CV(th=0.0001, ca=0.5, sd=1):
    """
    Finding the Coefficient Vector.

    Inputs:
    th -> threshold after which tail areas are considered to be insignificant.
    ca -> central coefficient. (By default equal to 0.5 because: "In fact, absence of
                                a specific choice of a0, a value of 0.5 for it would be
                                safe for most purposes, in the sense that the curve of
                                varying measure with time corresponding to a0 = 0.5 is
                                midway between the smoothest and the spikiest curves and
                                the average measure lies midway between the highest possible
                                and the lowest possible values. (see Fig 2,3)").
    sd -> standard deviation is = 1.

    Outputs:
    CV -> The coefficient vector.
    """

    gap = 2*norminv((ca + 1)/2, 0, sd)
    start = gap/2
    a = (ca + 1)/2
    b = ca

    while(1 - a) > th:
        c = normcdf(start + gap, 0, sd) - a
        a = normcdf(start + gap, 0, sd)
        b = horzcat(b, c)
        start = start + gap

    CV = horzcat(b[-1:0:-1], b)

    return CV



CV = find_CV(th=0.0001, ca=0.5, sd=1)


def calc_peak_sync_measure(pA, CV):
    """
    Calculating the measure of peak synchronization using the above found coefficient vector.

    Inputs:
    pA -> an r x N matrix, consisting of r many peak detected signals, with N time points each.
        r -> number of datasets with a (Signals)
    CV -> The coefficient vector.

    Ouputs:
    mc -> a 1 x N vector, with each index carrying the measure of peak synchronization for the
        r signals at that time point.
    """

    r = pA.shape[0]
    N = pA.shape[1]

    F = np.zeros_like(pA, dtype=float) # r x N matrix of zeroes. Important that values are float.
    n = int(np.floor(len(CV)/2))
    S = np.zeros((int(comb(r, 2)), N)) # rC2 matrix of zeroes

    for i in range(r):
        for t in range(n, N-n):
            F[i, t] = np.dot(pA[i, t-n:t+n+1], CV.transpose())
            # print(pA)
            # print(CV.transpose())

    cpt = 0 # z in paper. Starting at 0 here because of python indexes.
    for i in range(r):
        for j in range(i+1, r):
            for t in range(n, N-n):
                if np.dot(pA[i, t], pA[j,t]) == 0:
                    S[cpt, t] = max(np.dot(F[i,t], pA[j,t]), np.dot(F[j,t], pA[i,t]))
                else:
                    S[cpt, t] = (np.dot(F[i,t], pA[j,t]) + np.dot(F[j,t], pA[i,t]))/2
            cpt += 1

    mc = S.mean(axis=0)

    return mc








def peak_transform_sd(data, tm):
    """
    Using multiplier of standard-deviation.

    data: N x L matrix, N profiles of length L

    l = codon density signal profile sequence
    tm = threshold multipler.

    Returns individual profiles peak values (st) (Series object)
    """

    N, L = data.shape

    peak_mat = np.zeros(( N, L ))


    for i in range(N):
        current_data = data[i,:]
        
        for j in range(L):
            if current_data[j] >= np.mean(current_data) + tm*np.std(current_data) :
                peak_mat[i,j] = 1

    return peak_mat







def biswas_paper_peak_transform(data):
    """
    A. Peak Detection
    Transform peaks above treshold to 1. If not a peak, then 0.
    Method used by the Biswas et al. paper.

    Threshold = median + 2*sd

    window size = 2w+1  
    w = half window size
    """

    w = 7
    N, L = data.shape

    peak_mat = np.zeros(( N, L ))


    for i in range(N):

        current_data = data[i,:]
        current_mean = np.mean(current_data)
        
        for j in range(w, L-(w+1)):
            if current_data[j] >= current_mean + 2*np.std(current_data[(j-w):(j+w+1)]) :
                peak_mat[i,j] = 1

    return peak_mat




def run_mc(dataIN, CV):
    """
    wrapper
    """

    #current_peaks = biswas_paper_peak_transform(dataIN)
    current_peaks = peak_transform_sd(dataIN, tm=1)

    output_mc = calc_peak_sync_measure(current_peaks, CV)   

    return output_mc, current_peaks





def shuffle_data(orf_mat):
    """
    return matrix of randomized profiles
    """

    N = orf_mat.shape[0]

    rand_mat = np.copy(orf_mat)

    for i in range(N):
        local_state = np.random.RandomState(seed=None)
        rand_mat[i,:] = local_state.choice(rand_mat[i,:], len(rand_mat[i,:]), replace=False)
        #np.random.shuffle(rand_mat[i,:])

    return rand_mat




def get_rand_mc(dataIN):
    """
    get MC from randomized profile
    """
    current_matrix = np.copy(dataIN)
    current_matrix = shuffle_data(current_matrix)
    current_mc = run_mc(current_matrix, CV)

    return current_mc





def run_mc_frompeaks(peaksIN, CV=CV):
    """
    wrapper
    """

    #current_peaks = biswas_paper_peak_transform(dataIN)
    #current_peaks = peak_transform_sd(dataIN, tm=1)    

    output_mc = calc_peak_sync_measure(peaksIN, CV)   

    return output_mc



def rand_mc_frompeaks(peaksIN, CV=CV):
    """
    wrapper
    """

    current_peaks = np.copy(peaksIN)
    current_peaks2 = shuffle_data(current_peaks)
    current_mc = run_mc_frompeaks(current_peaks2, CV)
#    print(np.sum(current_peaks), np.sum(current_peaks2) )

    return current_mc




if __name__ == "__main__":

    
    with open("scikit_mat.pkl", 'rb') as f_ribo:
        scikit_data = pickle.load(f_ribo)

    with open("mm_consensus.pkl", 'rb') as f_mm:
        mm_consensus = pickle.load(f_mm)


    ### Compute consensus profile: ###
    # Compute corefficient vector (Probability Distribution):
    # th = 0.0001
    # ca = 0.5
    CV = find_CV(th=0.0001, ca=0.5, sd=1)

    list_orfs = list( scikit_data.keys() )

    mc_dict = {}
    theta_df = pd.DataFrame(columns=['ORF', 'p_5', 'p_10', 'p_20', 'p_80', 'p_90', 'p_95', 'p3_5', 'p3_10', 'p3_20', 'p3_80', 'p3_90', 'p3_95'])

    peaks = True

    counter = 0
    for ix, orf in enumerate(list_orfs):

        current_data = scikit_data[orf]

        current_mm = mm_consensus[orf]  # boolean: True for good sequence, False for multi-mapping
        print(ix, orf, current_data.shape[1], len(current_mm))
        if current_data.shape[1] == len(current_mm):


            current_data_mm = current_data[:,current_mm]        # for randomized consensus, chop         
            
            current_data[:,~current_mm] = 0                     # after, for false consensus (i.e. multimapping), set to 0
            mc_dict[orf], current_peaks = run_mc(current_data, CV)

            if peaks:
                max_iter = 100        
                pool = mp.Pool(processes=10)
                output = pool.map(rand_mc_frompeaks, [current_peaks for iteration in range(max_iter)])
                output = np.array(output)
                pool.close()

            else:
                max_iter = 100        
                pool = mp.Pool(processes=10)
                output = pool.map(get_rand_mc, [current_data_mm for iteration in range(max_iter)])
                output = np.array(output)
                pool.close()

        
            output3 = np.zeros(( output.shape[0], output.shape[1]-2 ))
            for rand_experiment in range(output3.shape[0]):
                for position in range(output3.shape[1]-2):     #to get kmers of length 3
                    output3[rand_experiment, position] = np.mean(output[rand_experiment, position:position+3])

            p_5  = np.around( np.percentile(output,  5), 5)
            p_10 = np.around( np.percentile(output, 10), 5)
            p_20 = np.around( np.percentile(output, 20), 5)
            p_80 = np.around( np.percentile(output, 80), 5)
            p_90 = np.around( np.percentile(output, 90), 5)
            p_95 = np.around( np.percentile(output, 95), 5)
  
            p3_5  = np.around( np.percentile(output3,  5), 5)
            p3_10 = np.around( np.percentile(output3, 10), 5)
            p3_20 = np.around( np.percentile(output3, 20), 5)
            p3_80 = np.around( np.percentile(output3, 80), 5)
            p3_90 = np.around( np.percentile(output3, 90), 5)
            p3_95 = np.around( np.percentile(output3, 95), 5)
  


            theta_df.loc[counter] = [orf, p_5, p_10, p_20, p_80, p_90, p_95, p3_5, p3_10, p3_20, p3_80, p3_90, p3_95]
            counter += 1


            # write intermediate output for debugging
            if counter % 10 == 0:
                theta_df.to_csv("theta.txt", header=True, index=False, sep='\t')
       
                with open('mc_dict.pkl', 'wb') as f:
                    pickle.dump(mc_dict, f, protocol=pickle.HIGHEST_PROTOCOL)



    with open('mc_dict.pkl', 'wb') as f:
        pickle.dump(mc_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    theta_df.to_csv("theta.txt", header=True, index=False, sep='\t')


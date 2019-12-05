import os, sys
import numpy as np
import pandas as pd
import subprocess
import glob
import json
import csv
import pickle
from scipy import stats
from Bio.Seq import Seq
from itertools import product
from peak_synchronization import *



#--------------------------------------------------------------------------------------------
def load_data():
    """
    load required base data
    """

    scikit_data = pickle.load( open("../data/processed/scikit_mat.pkl", 'rb') )

    PATH_MM_consensus = "../data/processed/mm_consensus.pkl"
    if os.path.exists(PATH_MM_consensus):
        mm_consensus = pickle.load( open(PATH_MM_consensus, 'rb') )
    else:
        mm_data = pickle.load( open("../data/processed/mm_mat.pkl", 'rb') )
        mm_consensus = mm_consensus(mm_data, PATH_MM_consensus)

    scikit_cons = pickle.load( open("../data/processed/mc_dict.pkl", 'rb') )

    return scikit_data, scikit_cons, mm_consensus



#--------------------------------------------------------------------------------------------
def mm_consensus(data_mm, fileOUT):
    """
    build consensus from multi-mapping data
    initial MM boolean is True for multi-mapping
    arbitrary threshold of MM in at least 5 indv. profiles for consesnsus MM
    converted to consensus MM boolean that is True for good position, False for MM position! 
    """

    mm_profile = {}
    theta_mm = 5
    for orf in data_mm.keys():
        current_mat = data_mm[orf]
        current_bool = np.sum(current_mat, 0) <= theta_mm
        mm_profile[orf] = current_bool

    pickle.dump(mm_profile, open(fileOUT))

    return mm_profile



#--------------------------------------------------------------------------------------------
def cor2consensus(data_scikit, data_consensus, data_mm):
    """
    compuare pairwise correlations between indv. profiles and consensus and their average
    """

    corr_df = pd.DataFrame(columns=['ORF', 'Correlation'])

    list_orfs = list( data_consensus.keys() )
    for ix, orf in enumerate(list_orfs):
        output = np.zeros(( 20 ))
        coef = 0
        p = 1

        current_data = data_scikit[orf]
        current_cons = data_consensus[orf]
        current_mm   = data_mm[orf]

        for i in range(20):
            data1 = current_data[i, current_mm]
            data2 = current_cons[current_mm]

            if np.sum(data1 > 0) > 10 and np.sum(data2 > 0) > 10: 		# at least 10 positions with signal
                coef, p = stats.spearmanr(current_data[i,current_mm], current_cons[current_mm])
                output[i] = coef
 
        corr = np.around( np.mean(output), 5)
        corr_df.loc[ix] = (orf, corr)

    return corr_df




#--------------------------------------------------------------------------------------------
def consensus_of_random(data):
    """
    assemble random sets of 'replicas':
    - for each gene, pick one random replica
    - prune by multi-mapping (only good parts)
    - get minimum length of profiles and assemble
    - compute consensus
    - compute avrg correlation to consensus
    """

    list_orfs = list( data.keys() )

    result = []

    max_iter = 2500  #about the same number as the set of genes analyzed 
    i = 0
    while i < max_iter:
        current_profiles = []
        current_profiles_Lmin = []
        current_lengths = []
      
        rand_genes = list( np.random.choice(list_orfs, 20) )

        for orf in rand_genes:
            current_rand = data[orf]
            current_profile = current_rand[ np.random.randint( 0,current_rand.shape[0] ) ]
            current_profiles.append(current_profile)
            current_lengths.append( len(current_profile) )

        Lmin = np.min( np.array(current_lengths) )
        
        if Lmin > 100:		# only generate random "replicas" of sufficient length, but doesn't really matter
            for j in range( len(current_profiles) ):
                current_profile = current_profiles[j]
                current_profile = current_profile[0:Lmin]
                current_profiles_Lmin.append(current_profile)
            current_profiles_Lmin = np.array(current_profiles_Lmin)

            CV = find_CV(th=0.0001, ca=0.5, sd=1)
            current_mc, current_peaks = run_mc(current_profiles_Lmin, CV)

            list_corr = []
            for j in range(20):
                data1 = current_profiles_Lmin[j,:]           
                data2 = np.copy(current_mc)
                N_data1 = np.sum(data1 > 0)
                N_data2 = np.sum(data2 > 0) 
                if (N_data1 > 10) and (N_data2 > 10): 		# at least 10 positions with signal, but no need to put actually
                    coef, p = stats.spearmanr(data1, data2)
                    list_corr.append(coef)
 
            corr = np.around( np.mean( np.array(list_corr)) , 5) 
            result.append(corr)
            i += 1

    return result





#--------------------------------------------------------------------------------------------
def data_figures1AB():
    """
    data for figure panels 1A and 1B
    """

    profiles_YKL143W = scikit_data['YKL143W']
    consensus_YKL143W = scikit_cons['YKL143W']

    np.savetxt("../data/figures/figure1/YKL143W_profiles.txt", profiles_YKL143W)
    np.savetxt("../data/figures/figure1/YKL143W_consensus.txt", consensus_YKL143W)


#--------------------------------------------------------------------------------------------
def data_figure1C():
    """
    data for figure panel 1C
    """

    corr_real = cor2consensus(scikit_data, scikit_cons, mm_consensus)
    corr_real.to_csv("../data/figures/figure1/corr2mc.txt", header=True, index=False, sep='\t')

    corr_rand = consensus_of_random(scikit_data)
    np.savetxt("../data/figures/figure1/corr2randmc.txt", np.array(corr_rand))


#--------------------------------------------------------------------------------------------
def data_figure1D():
    """
    data for figure panel 1C
    """

    list_orfs = list( scikit_cons.keys() )

    all_consensus = np.array([])

    for orf in list_orfs:
        current_consensus = scikit_cons[orf]
        current_mm = mm_consensus[orf]

        if len(current_consensus) == len(current_mm):

            current_consensus_mm = current_consensus[current_mm]
            current_consensus_mm = np.around(current_consensus, 4)
            all_consensus = np.concatenate((all_consensus, current_consensus_mm))

    np.savetxt("../data/figures/figure1/all_consensus.txt", all_consensus, fmt='%.4f')




#--------------------------------------------------------------------------------------------
def data_figure1E():
    """
    compare to expression data from processing/filtering step
    """

    yeast_filter = pd.read_csv("../data/processed/yeast_filter.txt", header=0, index_col=False, sep='\t')
    corr2mc = pd.read_csv("../data/figures/figure1/corr2mc.txt", header=0, index_col=False, sep='\t')

    peaks_df = pd.DataFrame(columns=['ORF', 'EL', 'corr', 'cons90'])

    list_orfs = list( scikit_cons.keys() )
    counter = 0
    for orf in list_orfs:
        current_data = scikit_data[orf]
        current_consensus = scikit_cons[orf]
        current_mm = mm_consensus[orf]

        current_el = yeast_filter[yeast_filter['ORF']==orf]['EL'].item()
        current_corr = corr2mc[corr2mc['ORF']==orf]['Correlation'].item()

        current_data = current_data[:, current_mm]
        current_consensus= current_consensus[current_mm] 

        th90_cons = np.around( np.percentile(current_consensus, 90), 5)

        peaks_df.loc[counter] = (orf, current_el, current_corr, th90_cons )
        counter += 1

    peaks_df.to_csv("../data/figures/figure1/data_qc.txt", header=True, index=False, sep='\t')







if __name__ == '__main__':


    scikit_data, scikit_cons, mm_consensus = load_data()

    data_figures1AB()
    data_figure1C()
    data_figure1D()
    data_figure1E()

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
import multiprocessing as mp
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
def cor2consensus_crossvalidated(data_scikit, data_mm):
    """
    leave-one-out crossvalidation of correlation to consenus
    """

    CV = find_CV(th=0.0001, ca=0.5, sd=1)
    corr_df = pd.DataFrame(columns=['ORF', 'Correlation'])
    list_orfs = list( data_scikit.keys() )
    for ix, orf in enumerate(list_orfs):
        output = np.zeros(( 20 ))
        coef = 0
        p = 1

        current_data = data_scikit[orf]
        current_mm   = data_mm[orf]         
    
        if np.shape(current_data)[1] == len(current_mm):

            current_data[:,~current_mm] = 0                     # after, for false consensus (i.e. multimapping), set to 0

            pool = mp.Pool(processes=10)
            output = pool.map(looxv_multi, [(current_data, current_mm, idx, CV) for idx in range(20)])
            output = np.array(output)
            pool.close()
 
            corr = np.around( np.nanmean(output), 5)
            corr_df.loc[len(corr_df)] = (orf, corr)

            print(ix, orf, corr)

    return corr_df



def looxv_multi(input_data):
    """
    subroutine to use multi-processing for leave-one-out crossvalidation
    """

    sub_data = input_data[0]
    sub_mm = input_data[1]
    sub_idx = input_data[2]
    CV = input_data[3]

    data1 = sub_data[sub_idx, sub_mm]
    sub_data_xv = sub_data[np.arange(np.shape(sub_data)[0])!=sub_idx, :]
    sub_cons_xv, sub_peaks = run_mc(sub_data_xv, CV)
    data2 = sub_cons_xv[sub_mm]

    if np.sum(data1 > 0) > 10 and np.sum(data2 > 0) > 10:       # at least 10 positions with signal
        coef, p = stats.spearmanr(data1, data2)
    else:
        coef = np.nan

    return coef



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
 
            corr = np.around( np.nanmean( np.array(list_corr)) , 5) 
            result.append(corr)
            i += 1

    return result




#--------------------------------------------------------------------------------------------
def consensus_of_random_xv(data):
    """
    assemble random sets of 'replicas':
    - for each gene, pick one random replica
    - prune by multi-mapping (only good parts)
    - get minimum length of profiles and assemble
    - compute consensus
    - compute avrg correlation to consensus, cross-validated by leave-one-out
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
        
        if Lmin > 100:      # only generate random "replicas" of sufficient length, but doesn't really matter
            for j in range( len(current_profiles) ):
                current_profile = current_profiles[j]
                current_profile = current_profile[0:Lmin]
                current_profiles_Lmin.append(current_profile)
            current_profiles_Lmin = np.array(current_profiles_Lmin)

            CV = find_CV(th=0.0001, ca=0.5, sd=1)
            #current_mc, current_peaks = run_mc(current_profiles_Lmin, CV)

            list_corr = []

            pool = mp.Pool(processes=10)
            output = pool.map(looxv_multi, [(current_profiles_Lmin, np.ones( (np.shape(current_profiles_Lmin)[1]), dtype=bool), idx, CV) for idx in range(20)])
            output = np.array(output)
            pool.close()
 
            corr = np.around( np.nanmean( np.array(output)) , 5) 
            result.append(corr)
            i += 1
    
            print(i, corr)

    return result



#--------------------------------------------------------------------------------------------
def peak_stats(data_scikit, data_cons, data_mm):
    """
    summary stats on consensus peaks
    """

    CV = find_CV(th=0.0001, ca=0.5, sd=1)
    theta_cons = 0.1                # arbitrary threshold to get some stats on peaks in consensus
    peak_df = pd.DataFrame(columns=['ORF', 'peaks_cons', 'peaks_mean', 'peaks_sd'])

    list_orfs = list( data_cons.keys() )
    for ix, orf in enumerate(list_orfs):
        output = np.zeros(( 20 ))
        coef = 0
        p = 1

        current_data = data_scikit[orf]
        current_mm   = data_mm[orf] 
        current_cons = data_cons[orf]        

        if np.shape(current_data)[1] == len(current_mm):
            current_data[:,~current_mm] = 0                     # after, for false consensus (i.e. multimapping), set to 0
            current_cons, current_peaks = run_mc(current_data, CV)
            current_num_peaks = np.sum(current_peaks, 1)
            current_peaks_cons = np.sum( current_cons > theta_cons )

            peak_df.loc[len(peak_df)] = (orf, current_peaks_cons, np.mean(current_num_peaks), np.around(np.std(current_num_peaks),3) )
            print(ix, orf, current_peaks_cons, np.mean(current_num_peaks), np.around(np.std(current_num_peaks),3) )

    return peak_df




#--------------------------------------------------------------------------------------------
def biswas_robustness(data_scikit, data_mm):
    """
    summary stats on consensus peaks
    """

    CV = find_CV(th=0.0001, ca=0.5, sd=1)
    CV_th001 = find_CV(th=0.001, ca=0.5, sd=1)
    CV_th01 = find_CV(th=0.01, ca=0.5, sd=1)
    CV_th00001 = find_CV(th=0.00001, ca=0.5, sd=1)
    CV_sd15 = find_CV(th=0.0001, ca=0.5, sd=1.5)
    CV_sd05 = find_CV(th=0.0001, ca=0.5, sd=0.5)
    CV_ca09 = find_CV(th=0.0001, ca=0.9, sd=0.5)
    CV_ca01 = find_CV(th=0.0001, ca=0.1, sd=0.5)


    biswas_df = pd.DataFrame(columns=['ORF', 'corr_th001', 'corr_th01', 'corr_th00001', 'corr_sd15', 'corr_sd05', 'corr_ca09', 'corr_ca01'])

    list_orfs = list( data_scikit.keys() )
    for ix, orf in enumerate(list_orfs):
        output = np.zeros(( 7 ))
        coef = 0
        p = 1

        current_data = data_scikit[orf]
        current_mm   = data_mm[orf]         

        if np.shape(current_data)[1] == len(current_mm):
            current_data[:,~current_mm] = 0                     # after, for false consensus (i.e. multimapping), set to 0

            current_cons, current_peaks = run_mc(current_data, CV)
            current_cons_th001, current_peaks_th001 = run_mc(current_data, CV_th001)
            current_cons_th01, current_peaks_th01 = run_mc(current_data, CV_th01)
            current_cons_th00001, current_peaks_th00001 = run_mc(current_data, CV_th00001)
            current_cons_sd15, current_peaks_sd15 = run_mc(current_data, CV_sd15)
            current_cons_sd05, current_peaks_sd05 = run_mc(current_data, CV_sd05)
            current_cons_ca09, current_peaks_ca09 = run_mc(current_data, CV_ca09)
            current_cons_ca01, current_peaks_ca01 = run_mc(current_data, CV_ca01)

            output[0], p = stats.spearmanr(current_cons, current_cons_th001)
            output[1], p = stats.spearmanr(current_cons, current_cons_th01)
            output[2], p = stats.spearmanr(current_cons, current_cons_th00001)
            output[3], p = stats.spearmanr(current_cons, current_cons_sd15)
            output[4], p = stats.spearmanr(current_cons, current_cons_sd05)
            output[5], p = stats.spearmanr(current_cons, current_cons_ca09)
            output[6], p = stats.spearmanr(current_cons, current_cons_ca01)
            output = np.around(output,3)

            biswas_df.loc[len(biswas_df)] = ( orf, output[0], output[1], output[2], output[3], output[4], output[5], output[6] ) 
            print(ix, orf, output[0], output[1], output[2], output[3], output[4], output[5], output[6] ) 

    return biswas_df





#--------------------------------------------------------------------------------------------
def parallel_mc(packaged_input):

    dataset = packaged_input[0]
    current_data = packaged_input[1]
    current_mm = packaged_input[2]
    current_cons = packaged_input[3]
    CV = packaged_input[4]


    coef = np.nan
    dataset_selector = np.ones( (np.shape(current_data)[0] ), dtype=bool)

    if len(np.shape(dataset)) == 0:
        dataset_selector[dataset] = False
    elif len(np.shape(dataset)) > 0:
        for i in list(dataset):
            dataset_selector[i] = False
    current_data = current_data[dataset_selector, :]   # leave out dataset

    if np.shape(current_data)[1] == len(current_mm):
        current_data[:,~current_mm] = 0     
        new_cons, new_peaks = run_mc(current_data, CV)
        if np.sum(current_cons[current_mm] > 0) > 10 and np.sum(new_cons[current_mm] > 0) > 10: 
            coef, p = stats.spearmanr(current_cons[current_mm], new_cons[current_mm])
  
    return coef



#--------------------------------------------------------------------------------------------
def dataset_robustness(data_scikit, data_cons, data_mm):
    """
    summary stats on consensus peaks
    """

    CV = find_CV(th=0.0001, ca=0.5, sd=1)
   
    list_orfs = list( data_cons.keys() )
    result_df = pd.DataFrame({'ORF': list_orfs})

    
    for dataset in range(20):
        print("processing dataset", dataset)
        correlation_result = np.zeros(( len(list_orfs) )) * np.nan

        pool = mp.Pool(processes=10)
        pool_input = [(dataset, data_scikit[list_orfs[i]], data_mm[list_orfs[i]], data_cons[list_orfs[i]], CV) for i in range(len(list_orfs))]
        output = pool.map(parallel_mc, pool_input)
        pool.close()

        correlation_result = np.around(np.array(output), 3)

        result_df[list_SRA[dataset]] = list(correlation_result)


    return result_df



#--------------------------------------------------------------------------------------------
def dataset_robustness_selected(data_scikit, data_cons, data_mm, datasets):
    """
    summary stats on consensus peaks
    """

    CV = find_CV(th=0.0001, ca=0.5, sd=1)
   
    list_orfs = list( data_cons.keys() )
    result_df = pd.DataFrame({'ORF': list_orfs})

    correlation_result = np.zeros(( len(list_orfs) )) * np.nan

    pool = mp.Pool(processes=10)
    pool_input = [(datasets, data_scikit[list_orfs[i]], data_mm[list_orfs[i]], data_cons[list_orfs[i]], CV) for i in range(len(list_orfs))]
    output = pool.map(parallel_mc, pool_input)
    pool.close()

    correlation_result = np.around(np.array(output), 3)

    result_df["Correlation"] = list(correlation_result)


    return result_df




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





#--------------------------------------------------------------------------------------------
def data_suppl_figure1():
    """
    data for supplementary figure panel
    """

#    print("computing cross-validated correlation to consensus")
#    corr_real_xv =  cor2consensus_crossvalidated(scikit_data, mm_consensus)
#    corr_real_xv.to_csv("../data/figures/figure1/suppl_corr2mc_xv.txt", header=True, index=False, sep='\t')

#    print("computing cross-validated correlation to random consensus")
#    corr_rand = consensus_of_random_xv(scikit_data)
#    np.savetxt("../data/figures/figure1/suppl_corr2randmc_xv.txt", np.array(corr_rand))

#    print("computing peak stats")
#    peak_df = peak_stats(scikit_data, scikit_cons, mm_consensus)
#    peak_df.to_csv("../data/figures/figure1/suppl_peak_stats.txt", header=True, index=False, sep='\t')

#    print("checking robustness of biswas consensus profile")
#    biswasdf = biswas_robustness(scikit_data, mm_consensus)
#    biswasdf.to_csv("../data/figures/figure1/suppl_biswas.txt", header=True, index=False, sep='\t')

    print("checking robustness of consensus with respect to indv. datasets")
    resultdf = dataset_robustness(scikit_data, scikit_cons, mm_consensus)
    resultdf.to_csv("../data/figures/figure1/suppl_datasets.txt", header=True, index=False, na_rep="NA", sep='\t')

    print("checking robustness without Lecanda data (1min of CHX)")
    lecandadf = dataset_robustness_selected(scikit_data, scikit_cons, mm_consensus, list([13,14]) )
    lecandadf.to_csv("../data/figures/figure1/suppl_lecanda.txt", header=True, na_rep="NA", index=False, sep='\t')





if __name__ == '__main__':


    scikit_data, scikit_cons, mm_consensus = load_data()

    list_SRA = ['SRR1002819','SRR1042855','SRR1520311','SRR1688545','SRR1944981','SRR1944982','SRR1944983',
                'SRR2046309','SRR2046310','SRR3493886','SRR3493887','SRR3623557','SRR3623558','SRR3945926',
                'SRR3945928','SRR4000288','SRR4000289','SRR5008134','SRR5008135','SRR5090936']


#    data_figures1AB()
#    data_figure1C()
#    data_figure1D()
#    data_figure1E()

    data_suppl_figure1()
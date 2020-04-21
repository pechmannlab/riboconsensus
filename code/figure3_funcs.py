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
def compute_theta():
    """
    get thresholds from randomized consensus
    """

    ### Compute consensus profile: ###
    # Compute corefficient vector (Probability Distribution):
    # th = 0.0001
    # ca = 0.5
    CV = find_CV(th=0.0001, ca=0.5, sd=1)

    list_orfs = list( scikit_data.keys() )
    mc_dict = {}
    theta_df = pd.DataFrame(columns=['ORF', 'p_5', 'p_10', 'p_20', 'p_80', 'p_90', 'p_95', 'p3_5', 'p3_10', 'p3_20', 'p3_80', 'p3_90', 'p3_95'])

    peaks = True        # reshuffle previously identified peak positions instead of consernsus profile as here equivalent and faster

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

    theta_df.to_csv("../data/figures/figure3/theta.txt", header=True, index=False, sep='\t')



#--------------------------------------------------------------------------------------------
def count_kmer(gene_list, codon_seqs, R, kmer_size=3):
    """
    loop through yeast transcriptome to get general expected redundancies of codon kmers
    this will generate a smallish list of candidate kmers within seconds, so speeds up search

    codon_seqs: dictionary of codon sequence as array of codons
    R: number of occurrances in genome
    kmer_size: length of kmer in number of codons

    MM: 'yes', if multi-mapping, exclude mm regions and count; 'no' count all
    """

    kmer = kmer_size
    MM = 'yes'

    list_seqfile = list( codon_seqs.keys() )
    kmer_dict = {}

    for orf in gene_list:
        if orf in list_seqfile:
            current_seq = np.array(codon_seqs[orf])

            for pos in range(len(current_seq) - (kmer + 1) ):
                if MM == 'yes' and orf in list( mm_consensus.keys() ):
                    current_mm = mm_consensus[orf]
                    if np.all(current_mm[pos:(pos+kmer)]):              # check that no kmer position is MM
                        current_kmer = "".join( current_seq[pos:pos+kmer])
                        if current_kmer in kmer_dict.keys():
                            kmer_dict[current_kmer] += 1
                        else:
                            kmer_dict[current_kmer] = 1

                elif MM == 'no':
                    current_kmer = "".join( current_seq[pos:pos+kmer])
                    if current_kmer in kmer_dict.keys():
                        kmer_dict[current_kmer] += 1
                    else:
                        kmer_dict[current_kmer] = 1

    new_dict = {}
    list_redundant = []
    for k in kmer_dict.keys():
        if kmer_dict[k] > R:
            if k not in list_redundant:
        	    list_redundant.append(k)
  
    return list_redundant

 




#--------------------------------------------------------------------------------------------
def extract_kmers(data_mc, data_mm, codon_seqs):
    """
    extracts all kmers of length 3 codons as per following criteria:
    - redundancy of at least 20 occurrences in set of analyzed genes
    - using pre-computed per-gene theta-thresholds ('theta.txt')
    - at least 10 uccurrences below and 10 above gene-specific thresholds
    """

    theta = pd.read_csv("../data/figures/figure3/theta.txt", header=0, index_col=False, sep='\t')

    trim = 20			# omit first and last 20 codons per gene due to known biases of translation initiation and termination
    kmer = 3

    list_orfs = list( data_mc.keys() )
    list_candidates = count_kmer(list_orfs, codon_seqs, 20, kmer)

    kmer_df = pd.DataFrame(columns=['ORF', 'kmer', 'position', 'class', 'threshold', 'all'])

    for ix, orf in enumerate(list_orfs):
        print(ix, orf)
        current_consensus = data_mc[orf]
        current_sequence = codon_seqs[orf]
        current_mm = data_mm[orf]

        current_mean = np.mean(current_consensus[current_mm])
        current_std  = np.std(current_consensus[current_mm])
        current_bottom5 = np.percentile( current_consensus[current_mm], 5 ) 

        # change significance thresholds here if wanting to test different ones! 
        current_theta_lo =  theta[theta['ORF']==orf]['p3_10'].item()
        current_theta_hi =  theta[theta['ORF']==orf]['p3_90'].item()    

        for pos in range( trim, len(current_sequence) - (trim + kmer) ):	# omit first/last 20 positions and allow for kmer length
            current_kmer = current_sequence[pos:pos+kmer]
            current_kmer2 = "".join(current_kmer)

            if current_kmer2 in list_candidates and np.all(current_mm[pos:pos+kmer]):								# omit multi-mapping positions
                current_score = current_consensus[pos:pos+kmer]
                current_threshold = 0
                current_all = 1
                current_class = np.nan

                if np.mean( current_score) <= current_theta_lo:
                    current_class = 0
                    current_threshold = 1
                  
                elif np.mean( current_score) > current_theta_hi:
                    current_class = 1
                    current_threshold = 1
                    
                else:
                    current_class = np.nan

                kmer_df.loc[len(kmer_df)] = (orf, current_kmer2, pos, current_class, current_threshold, current_all)
    kmer_df.to_csv("../data/figures/figure3/kmer_all.txt", header=True, index=False, sep='\t')


    # filter hits
    def filter_kmertable(kmertable, theta_min=10):
        """
        subfunction to only keep kmers with 10 cases of each lo/hi
        theta_min = 10 each (below / above)
        """
        kmertable_filtered = kmertable[kmertable['threshold']==1]

        list_kmers = list(set(kmertable_filtered['kmer']))
        list_keep = np.array([])

        for kmer in list_kmers:
            current_df = kmertable_filtered[kmertable_filtered['kmer']==kmer]
            if len(current_df) >= (2.*theta_min) :
                current_lo = current_df[current_df['class'] == 0]
                current_hi = current_df[current_df['class'] == 1]
                if len(current_lo) >= theta_min and len(current_hi) >= theta_min:
                    current_idx = np.where( np.array( kmertable_filtered['kmer']==kmer) )[0]
                    list_keep = np.concatenate((list_keep, current_idx))

        Fout_result = '../data/figures/figure3/kmer_filtered.txt'
        result_df = kmertable_filtered.iloc[list_keep]
        result_df.to_csv(Fout_result, header=True, index=False, sep='\t')

        return result_df

    df_filtered = filter_kmertable(kmer_df, 10)


    def filter_kmertable_nonDT(kmertable):
        """
        subfunction to only keep kmers with 15/5 cases of each lo/hi
        col_id: 'theta0595' or 'theta1090'
        hi: > 16 hi and < 4 lo
        lo: > 16 lo and < 4 hi
        """

        theta_min = 10
        theta_hi = 16
        theta_lo = 4

        # only kmers below/above gene-specific thresholds
        kmertable_filtered = kmertable[kmertable['threshold']==1]

        list_kmers = list(set(kmertable_filtered['kmer']))
        list_keep_hi = np.array([])
        list_keep_lo = np.array([])

        for kmer in list_kmers:
            current_df = kmertable_filtered[kmertable_filtered['kmer']==kmer]
            if len(current_df) >= (2.*theta_min) :
                current_lo = current_df[current_df['class'] == 0]
                current_hi = current_df[current_df['class'] == 1]

                if len(current_hi) >= theta_hi and len(current_lo) <= theta_lo:
                    current_idx = np.where( np.array( kmertable_filtered['kmer']==kmer) )[0]
                    list_keep_hi = np.concatenate((list_keep_hi, current_idx))

                elif len(current_hi) <= theta_lo and len(current_lo) >= theta_hi:
                    current_idx = np.where( np.array( kmertable_filtered['kmer']==kmer) )[0]
                    list_keep_lo = np.concatenate((list_keep_lo, current_idx))

        list_keep = np.concatenate( (list_keep_hi, list_keep_lo) )

        Fout_result = '../data/figures/figure3/kmer_filtered_nonDT+.txt'
        result_df_nonDT_hi = kmertable_filtered.iloc[list_keep_hi]
        result_df_nonDT_hi.to_csv(Fout_result, header=True, index=False, sep='\t')

        Fout_result = '../data/figures/figure3/kmer_filtered_nonDT-.txt'
        result_df_nonDT_lo = kmertable_filtered.iloc[list_keep_lo]
        result_df_nonDT_lo.to_csv(Fout_result, header=True, index=False, sep='\t')

        return result_df_nonDT_hi, result_df_nonDT_lo 

    df_filtered_nonDT_hi, df_filtered_nonDT_lo = filter_kmertable_nonDT(kmer_df)


    # compile count table
    kmer_counts = pd.DataFrame(columns=['kmer', 'lo', 'hi', 'tot', 'class'])
    for kmer in list(set(kmer_df['kmer'])):
        current_lo = 0
        current_hi = 0
        current_class = 'none'

        if kmer in list(set(df_filtered['kmer'])) :
        	current_class = 'DT' 
        	current_df = df_filtered[df_filtered['kmer']==kmer] 


        elif kmer in list(set(df_filtered_nonDT_hi['kmer'])):
            current_class = "nonDT+"
            current_df = df_filtered_nonDT_hi[df_filtered_nonDT_hi['kmer']==kmer]
        elif kmer in list(set(df_filtered_nonDT_lo['kmer'])):
            current_class = "nonDT-"
            current_df = df_filtered_nonDT_lo[df_filtered_nonDT_lo['kmer']==kmer]
        else:
            current_class = 'all'
            current_df = kmer_df[kmer_df['kmer']==kmer]

        current_lo = len( current_df[current_df['class'] == 0] )
        current_hi = len( current_df[current_df['class'] == 1] )
    
        kmer_counts.loc[len(kmer_counts)] = (kmer, current_lo, current_hi, current_lo + current_hi, current_class)
        kmer_counts.to_csv('../data/figures/figure3/kmer_counts.txt', header=True, index=False, sep='\t')
  
    return kmer_df, df_filtered, df_filtered_nonDT_hi, df_filtered_nonDT_lo




#--------------------------------------------------------------------------------------------
def kmer_frequencies(kmertable_all, kmertable_filtered, kmertable_nonDT_hi, kmertable_nonDT_lo, data_mm, codon_seqs):
    """
    codon frequencies in kmers
    background: mRNA sequences of set of 'good' genes (minus multi-mapping positions)
    redundant: all kmers with 20x redundancy in 'good' genes, i.e. 'all' from extract kmer function
    kmer_filtered: hits at 10%/90% thresholds
    """

    def codon_bgfreq(codon_seqs, data_mm):
        """
        get codon background frequencies from mRNA seqs
        seqs: dictionary of yeast mRNA sequences
        data_mc: dictionary of multi-mapping boolean
        """
        codon_counts = np.zeros(( len(codons_nonstop) ))
        list_orfs = list( data_mm.keys() )

        for ix, orf in enumerate(list_orfs):
            current_seq = codon_seqs[orf]
            current_mm = data_mm[orf]

            for pos in range( len(current_mm) ):
                if current_mm[pos] and current_seq[pos] in codons_nonstop:
                    current_index = codons_nonstop.index(current_seq[pos])
                    codon_counts[current_index] += 1
        codon_counts = np.around( codon_counts / np.sum(codon_counts), 5)

        return codon_counts


    def codonfreqs_kmerdf(kmertable):
        """
        get codon frequencies from kmertable
        """      
        codon_counts_kmer = np.zeros(( len(codons_nonstop) ))
        for kmer in kmertable['kmer']:
            current_kmer_codons = [ kmer[(i*3):((i*3)+3)] for i in range(3) ] # ! hard coded for length L=3
            for codon in current_kmer_codons:
                current_index = codons_nonstop.index(codon)
                codon_counts_kmer[current_index] += 1        
        codon_counts_kmer /= np.sum(codon_counts_kmer)

        return np.around(codon_counts_kmer, 5)

    #kmertable_threshold = kmertable_all[kmertable_all['threshold']==1]
    kmertable_all2      = kmertable_all[kmertable_all['threshold']==0]


    cc_bg = codon_bgfreq(codon_seqs, data_mm)
    cc_all  = codonfreqs_kmerdf(kmertable_all2)			# without hits
    cc_theta = codonfreqs_kmerdf(kmertable_filtered)
    cc_nDT_hi = codonfreqs_kmerdf(kmertable_nonDT_hi)   # min 16 max 4 at 1090
    cc_nDT_lo = codonfreqs_kmerdf(kmertable_nonDT_lo)   # min 16 max 4 at 1090

    output = pd.DataFrame({'codon': list(codons_nonstop), 
                            'kmer_theta': list(cc_theta), 
                            'redundant': list(cc_all), 
                            'background': list(cc_bg),
                            'nDThi': list(cc_nDT_hi),
                            'nDTlo': list(cc_nDT_lo) } )  
    output.to_csv("../data/figures/figure3/kmer_frequencies.txt", header=True, index=False, sep='\t')

    return output

  
#--------------------------------------------------------------------------------------------
def kmer_redundancy(list_orfs, yeast_codons):

    params = np.arange(3,8)
    result5 = np.zeros(( len(params) ))
    result10 = np.zeros(( len(params) ))
    result20 = np.zeros(( len(params) ))
    for ix, i in enumerate( params ):
        current_kmer5 = count_kmer( list_orfs, yeast_codons, 5, i)
        result5[ix] = len(current_kmer5)
        current_kmer10 = count_kmer( list_orfs, yeast_codons, 10, i)
        result10[ix] = len(current_kmer10)
        current_kmer20 = count_kmer( list_orfs, yeast_codons, 20, i)
        result20[ix] = len(current_kmer20)

    result_df = pd.DataFrame({"length": list(params), "counts5": list(result5), "counts10": list(result10), "counts20": list(result20)})
    result_df.to_csv("../data/figures/figure3/kmer_redundancy.txt", header=True, index=False, sep='\t')
    print(result_df)





if __name__ == '__main__':

 
    # load base data
    scikit_data = pickle.load(open("../data/processed/scikit_mat.pkl", 'rb'))
    scikit_consensus = pickle.load(open("../data/processed/mc_dict.pkl", 'rb'))
    mm_consensus = pickle.load(open("../data/processed/mm_consensus.pkl", 'rb'))
    sequence_codons = pickle.load(open('../data/processed/yeast_codons.pkl', 'rb'))


    # define codon list as global var
    codons_nonstop = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
    		  'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    		  'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    		  'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',]
   
 

    # ~~~~~~~ FUNCTIONS START HERE

    compute_theta()

    kmer_redundancy(list(scikit_data.keys()), sequence_codons)

    kmer_all, kmer_DT, kmer_nonDT_hi, kmer_nonDT_lo = extract_kmers(scikit_consensus, mm_consensus,  sequence_codons)

    kmer_freqs = kmer_frequencies(kmer_all, kmer_DT, kmer_nonDT_hi, kmer_nonDT_lo, mm_consensus, sequence_codons)

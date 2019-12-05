import os, sys
import numpy as np
import pandas as pd
import subprocess
import glob
import json
import csv
import pickle
from Bio.Seq import Seq
from itertools import product



#--------------------------------------------------------------------------------------------
def codon_translationrates_indprofiles(data_scikit, data_mm, codon_seq):
    """
    extracts distributions of codon translation rates from scikit ribo data
    scikit_data: dictionary of scikit_data per gene with profiles of each experiment
    codon_sequence: dictionary of mRNA sequence in codons
    """

    trim = 20 	# omit first and last 20 codons

    codon_duplets = [(x+y) for x in codons_nonstop for y in codons_nonstop]
  
    ALL_TR = []
    ALL2_TR = []
    TR_df = pd.DataFrame(columns=codons_nonstop)
    TR2_df = pd.DataFrame(columns=codon_duplets)
    TR2_raw_df = pd.DataFrame(columns=['codonpair', 'class', 'RD'])
    codon_counts = np.zeros(( 20, len(codons_nonstop) ))
    codon_duplet_counts = np.zeros(( 20, len(codon_duplets) ))

    list_orfs = list( data_scikit.keys() )

    counter_cp = 0

    for experiment in range( 20 ):
        codon_TR = [[] for i in range( len(codons_nonstop) )]
        codon2_TR = [[] for i in range( len(codon_duplets) )]
        print("Analyzing experiment", experiment)      
        for ix, orf in enumerate(list_orfs):
            current_data = data_scikit[orf]
            current_mm = data_mm[orf]
            current_seq = np.array( codon_seq[orf] )

            if current_data.shape[1] == len(current_mm):
                for pos in range(trim,  len(current_seq) - (trim+2) ):                             
                    if current_mm[pos]:  
                        current_codon = current_seq[pos]
                        current_codon_idx = codons_nonstop.index(current_codon)
                        current_TR = current_data[experiment,pos]

                        codon_TR[current_codon_idx].append(current_TR)   
                        codon_counts[experiment, current_codon_idx] += 1
                       
                    if current_mm[pos] and current_mm[pos+1]:
                        current_codon_pair = current_seq[pos] + current_seq[pos+1]
                        current_codon_pair_idx = codon_duplets.index(current_codon_pair)
                        current_pair_TR = (float(current_data[experiment, pos]) + float(current_data[experiment, pos+1])  )/2.

                        codon2_TR[current_codon_pair_idx].append(current_pair_TR)
                        codon_duplet_counts[experiment, current_codon_pair_idx] += 1

                        if current_codon_pair in good_inhibitory:
                            TR2_raw_df.loc[len(TR2_raw_df)] = (current_codon_pair, 'ind', current_pair_TR)
                        else:
                            counter_cp += 1
                            if counter_cp % 10 == 0:		# thin out 1 in 10 to reduce file size -> thin out more for faster run time!
                                TR2_raw_df.loc[len(TR2_raw_df)] = ('non', 'ind', current_pair_TR)
           
        TR_mean = [ np.around( np.mean(np.array(codon_TR[x])), 5)  for x in range( len(codons_nonstop) ) ]
        TR_median = [ np.around( np.median(np.array(codon_TR[x])), 5)  for x in range( len(codons_nonstop) ) ]
        TR_df.loc[experiment] = TR_mean

        TR2_mean = [ np.around( np.mean(np.array(codon2_TR[x])), 3)  for x in range( len(codon_duplets) ) ]
        TR2_df.loc[experiment] = TR2_mean
        
    TR_df.to_csv("../data/figures/figure2/codon_rates.txt", header=True, index=False, sep='\t')  
    TR2_df.to_csv("../data/figures/figure2/codon_duplets_rates.txt", header=True, index=False, sep='\t')

    TR2_raw_df.to_csv("../data/figures/figure2/codonpair_rates_raw10_ind.txt", header=True, index=False, sep='\t')

    np.savetxt("../data/figures/figure2/codon_counts.txt", codon_counts, fmt='%i')
    np.savetxt("../data/figures/figure2/codon_duplet_counts.txt", codon_duplet_counts, fmt='%i')

    return TR_df


#--------------------------------------------------------------------------------------------
def codon_translationrates_consensus(data_consensus, data_mean, data_mm, codon_seq):
    """
    extracts distributions of codon translation rates from scikit ribo data
    scikit_data: dictionary of scikit_data per gene with profiles of each experiment
    codon_sequence: dictionary of mRNA sequence in codons
    """

    trim = 20   # omit first and last 20 codons


    list_orfs = list( data_consensus.keys() )
    codon_TR = [[] for i in range( len(codons_nonstop) )]
    codon_TR_naive = [[] for i in range( len(codons_nonstop) )]

    codon_duplets = [(x+y) for x in codons_nonstop for y in codons_nonstop]
    codon2_TR = [[] for i in range( len(codon_duplets) )]
    codon2_TR_naive = [[] for i in range( len(codon_duplets) )]

    TR2_cons_df = pd.DataFrame(columns=['codonpair', 'class', 'RD'])

    counter_cp = 0

    for ix, orf in enumerate(list_orfs):
        current_consensus = data_consensus[orf]
        current_naivemean = data_mean[orf]
        current_mm = data_mm[orf]
        current_seq = np.array( codon_seq[orf] )
        print(ix, orf)
        if len(current_consensus) == len(current_mm):
            for pos in range(trim,  len(current_seq) - (trim+2) ):      # +2 so that can also extract codon pairs
                                  
                if current_mm[pos]:  
                    current_codon = current_seq[pos]
                    current_codon_idx = codons_nonstop.index(current_codon)
                    current_TR = np.around( current_consensus[pos], 5)
                    current_TR_naive = np.around( current_naivemean[pos], 5)
                    codon_TR[current_codon_idx].append(current_TR) 
                    codon_TR_naive[current_codon_idx].append(current_TR_naive)  

                if current_mm[pos] and current_mm[pos+1]:
                    current_codon_pair = current_seq[pos] + current_seq[pos+1]
                    current_codon_pair_idx = codon_duplets.index(current_codon_pair)
                    current_pair_TR = np.around( (float(current_consensus[pos]) + float(current_consensus[pos+1])  )/2., 5)
                    current_pair_TR_naive = np.around( (current_naivemean[pos] + current_naivemean[pos+1])/2., 5 )
                    codon2_TR[current_codon_pair_idx].append(current_pair_TR)
                    codon2_TR_naive[current_codon_pair_idx].append(current_pair_TR_naive)


                    if current_codon_pair in good_inhibitory:
                        TR2_cons_df.loc[len(TR2_cons_df)] = (current_codon_pair, 'mc', current_pair_TR)
                        TR2_cons_df.loc[len(TR2_cons_df)] = (current_codon_pair, 'mn', current_pair_TR_naive)
                    else:
                        counter_cp += 1
                        if counter_cp % 10 == 0:	# thin out 1 in 10 to reduce file size
                            TR2_cons_df.loc[len(TR2_cons_df)] = ('non', 'mc', current_pair_TR)
                            TR2_cons_df.loc[len(TR2_cons_df)] = ('non', 'mn', current_pair_TR_naive)

         
    TR2_cons_df.to_csv("../data/figures/figure2/codonpair_rates_raw10_cons.txt", header=True, index=False, sep='\t')


    codon_df = pd.DataFrame(columns=['Codon', 'p_25_mc', 'median_mc', 'p_75_mc', 'p_25_mn', 'median_mn', 'p_75_mn', 'tAI'])
    for ix, codon in enumerate(codons_nonstop):
        mc_median = np.around( np.percentile(codon_TR[ix], 50), 5)
        mc_p25 = np.around( np.percentile(codon_TR[ix], 25), 5)
        mc_p75 = np.around( np.percentile(codon_TR[ix], 75), 5) 

        mn_median = np.around( np.percentile(codon_TR_naive[ix], 50), 5)
        mn_p25 = np.around( np.percentile(codon_TR_naive[ix], 25), 5)
        mn_p75 = np.around( np.percentile(codon_TR_naive[ix], 75), 5) 

        current_tAI = tAI[ix]

        codon_df.loc[ix] = (codon, mc_p25, mc_median, mc_p75, mn_p25, mn_median, mn_p75, current_tAI)
                   
    codon_df.to_csv("../data/figures/figure2/codon_rates_consensus.txt", header=True, index=False, sep='\t')


    codon_pair_df = pd.DataFrame(columns=['Codonpair', 'p25', 'median', 'p75', 'mn_p25', 'mn_median', 'mn_p75',])
    for ix, codonpair in enumerate(codon_duplets):

        if len(codon2_TR[ix]) > 10:
            current_median = np.around( np.percentile(np.array(codon2_TR[ix]), 50), 5)
            current_p25 = np.around( np.percentile(np.array(codon2_TR[ix]), 25), 5)
            current_p75 = np.around( np.percentile(np.array(codon2_TR[ix]), 75), 5)

            current_mn_median = np.around( np.percentile(np.array(codon2_TR_naive[ix]), 50), 5)
            current_mn_p25 = np.around( np.percentile(np.array(codon2_TR_naive[ix]), 25), 5)
            current_mn_p75 = np.around( np.percentile(np.array(codon2_TR_naive[ix]), 75), 5)

        else:
            current_median = np.nan
            current_p25 = np.nan
            current_p75 = np.nan
            current_mn_median = np.nan
            current_mn_p25 = np.nan
            current_mn_p75 = np.nan

        codon_pair_df.loc[ix] = (codonpair, current_p25, current_median, current_p75, current_mn_p25, current_mn_median, current_mn_p75)

    codon_pair_df.to_csv("../data/figures/figure2/codonpair_rates_consensus.txt", header=True, index=False, sep='\t')


    return codon_df




def data_figure2():
    """
    generate data for figure 2 panels
    """
 
    individual_df = codon_translationrates_indprofiles(scikit_data, mm_consensus, sequence_codons)
    consensus_df = codon_translationrates_consensus(scikit_consensus, scikit_mean, mm_consensus, sequence_codons)

    individual_dfT = individual_df.T
    final_df = consensus_df[['Codon', 'median_mc', 'median_mn']]
    for i in range(20):
        current_label = 'R'+str(i+1)
        current_vals = np.around(np.array(individual_dfT[i]),5) 
        final_df.insert(final_df.shape[1], current_label, current_vals)
    final_df.insert(final_df.shape[1], 'tAI', tAI)
  
    final_df.to_csv("../data/figures/figure2/codonrates_all.txt", header=True, index=False, sep='\t')





if __name__ == '__main__':


    codons_nonstop = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
              'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
              'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
              'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',]

    tAI = [0.431, 0.616, 1, 0.27, 0.246, 0.488, 0.14, 0.677, 0.677, 0.123, 0.278, 0.054, 0.123, 0.576, 0.616, 0.8, 0.554, 0.431, 0.239, 0.189, 0.616, 0.089, 0.197, 0.123, 0, 0.266, 0.062, 0.369, 0.185, 0.062, 0.059, 0.027, 0.862, 0.985, 0.399, 0.433, 0.308, 0.488, 0.099, 0.677, 0.185, 0.985, 0.182, 0.433, 0.123, 0.621, 0.163, 0.862, 0.493, 0.216, 0.185, 0.488, 0.121, 0.677, 0.246, 0.369, 0.108, 0.431, 0.616, 0.754, 0.27]

    inhibitory = ['CGACCG', 'CGAGCG', 'CGAATA', 'CTCCCG', 'CGACGA', 'CGACGG', 'CGACTG', 'GTACCG', 'GTGCGA', 'GTACGA', 'CTGCCG', 'CTGCGA', 'ATACGA', 'ATACGG', 'AGGCGG', 'CTGATA', 'AGGCGA']
    good_inhibitory = [ "CGAATA", "GTACCG", "GTGCGA", "GTACGA", "CTGCCG", "CTGCGA", "ATACGA", "ATACGG", "AGGCGG", "CTGATA", "AGGCGA"]
      
    scikit_data = pickle.load(open("../data/processed/scikit_mat.pkl", 'rb'))
    mm_consensus = pickle.load(open("../data/processed/mm_consensus.pkl", 'rb'))       
    scikit_consensus = pickle.load(open("../data/processed/mc_dict.pkl", 'rb'))
    scikit_mean = pickle.load(open("../data/processed/scikit_mean.pkl", 'rb'))
    sequence_codons = pickle.load(open('../data/processed/yeast_codons.pkl', 'rb'))


    data_figure2()

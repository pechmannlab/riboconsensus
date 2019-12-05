import os, sys
import numpy as np
import pandas as pd
import pickle
from gprofiler import GProfiler



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def split_kmerdf(kmertable, codonlist):
    """
    split kmer DF into kmers with and without mod-codon
    hard-coded for kmers length 3codons/9nts but easy to change
    """

    modcodon = codonlist #['GAA'] 		# codon matching tRNA with U34 modification based on Fig. 3D
    theta_modcod = 1

    n = 3 									# length 3 codons

    list_with = []
    list_without = []

    for ix, kmer in enumerate(kmertable['kmer']):
        current_codons = [ kmer[(i*3):((i*3)+3)] for i in range(n) ]
        current_modcods = np.zeros(( 3 ))
        for jx, codon in enumerate(current_codons):
            if codon in modcodon:
                current_modcods[jx] = 1
        if np.sum(current_modcods) >= theta_modcod:
            list_with.append(ix)
        elif np.sum(current_modcods) == 0: #else:
            list_without.append(ix)

    kmerdf_with = kmertable.iloc[list_with]
    kmerdf_without = kmertable.iloc[list_without]

    return kmerdf_with, kmerdf_without




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def count_kmer(codon_seqs, data_mm, R=20, kmer_size=3):
    """
    TO DO: import this function from Figure 3 file. 
    loop through yeast transcriptome to get general expected redundancies of codon kmers
    this will generate a smallish list of candidate kmers within seconds, so speeds up search
    codon_seqs: dictionary of codon sequence as array of codons
    R: number of occurrances in genome
    kmer_size: length of kmer in number of codons
    MM: 'yes', if multi-mapping, exclude mm regions and count, 'no' count all
    """

    kmer = kmer_size
    MM = 'yes'

    list_seqfile = list( codon_seqs.keys() )

    kmer_dict = {}

    gene_list = list(data_mm.keys())
    for orf in gene_list:
        if orf in list_seqfile:
            current_seq = np.array(codon_seqs[orf])
            for pos in range(len(current_seq) - (kmer + 1) ):
                if MM == 'yes' and orf in list( mm_consensus.keys() ):
                    current_mm = mm_consensus[orf]
                    if np.all(current_mm[pos:(pos+kmer)]):
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def count_redundant(genelist, codon_seqs, list_redundant):
    """
    loop through codon sequence for each gene in genelist
    count number of positions / substrings that are in list_redundant
    """

    result = {}

    trim = 20

    for ix, orf in enumerate(genelist):
        #print(ix, orf)        
        current_seq = codon_seqs[orf]
        current_count = 0
        for pos in range(trim, len(current_seq) - trim - 2 ):
            current_kmer = "".join(current_seq[pos:(pos+3)])
            if current_kmer in list_redundant:
                current_count += 1
        result[orf] = current_count

    return result



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def kmer_density(kmertable, codon_seqs, data_mm):
    """
    get counts/density of DT sequences per gene
    """

    list_orf_all = list( data_mm.keys() )
    list_orf_kmer = list(set(kmertable['ORF']))

    list_redundant = count_kmer(codon_seqs, data_mm, R=20, kmer_size=3)
    all_redundant = count_redundant(list_orf_all, codon_seqs, list_redundant)

    counts = pd.DataFrame(columns=['ORF', 'redundant', 'total', 'high', 'low', 'length'])
    
    for orf in list_orf_all:

        current_total = 0
        current_high = 0
        current_low = 0

        current_length = len(codon_seqs[orf])
        current_redundant = all_redundant[orf]
        current_table = kmertable[kmertable['ORF']==orf]
        current_total = len( current_table )
        current_high = len( current_table[current_table['class']==1] )
        current_low  = len( current_table[current_table['class']==0] )
        
        counts.loc[len(counts)] = (orf, current_redundant, current_total, current_high, current_low, current_length ) 

    #counts.sort_values(['total'], ascending=False, inplace=True)


    return counts



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def get_genelist(counttable, col1, col2, method="difference"):
    """
    get list of genes from tails of distribution
    - col1 normalized (divided) by col2
    - loop through to deal with avoiding division by zero, inf, etc. 
    - 5% tails or top/bottom 100, whatever is more
    method: "difference" or "division"
    """

    score = np.zeros(( len(counttable) ))

    for i in range( len(counttable) ):
        current_orf = counttable.iloc[i]['ORF']
        current_col1 = int( counttable.iloc[i][col1] )
        current_col2 = int( counttable.iloc[i][col2] )

        if method == "difference":
            if ((current_col1 + current_col2)/2.) == 0:
                current_score = 0
            else:
                current_score = (current_col1 - current_col2) / ((current_col1 + current_col2)/2.)
            score[i] = current_score

        elif method == "division":
            if current_col1 == 0:
                score[i] = np.nan
            elif current_col2 == 0:
                score[i] = np.nan
            else:
                score[i] = current_col1 / current_col2

    counttable['score'] = score

    sel_notnan = np.isnan(score) == False  
    resultDF = counttable.loc[sel_notnan]
    resultDF = resultDF[['ORF', 'score']]

    # gets lists by thresholds and extremes
    theta_hi = np.percentile(resultDF['score'], 95)
    theta_lo = np.percentile(resultDF['score'], 5 )
    list_hi_theta = list( resultDF[ resultDF['score'] >= theta_hi ]['ORF'] )
    list_lo_theta = list( resultDF[ resultDF['score'] <= theta_lo ]['ORF'] )

    resultDF_sorted = resultDF.sort_values('score', ascending=False, inplace=False)
    L = len(resultDF_sorted)
    list_hi_ex = resultDF_sorted['ORF'][0:100]
    list_lo_ex = resultDF_sorted['ORF'][(L-100):L]

    if len(list_hi_theta) > len(list_hi_ex):
        list_hi = list_hi_theta
    else:
        list_hi = list_hi_ex

    if len(list_lo_theta) > len(list_lo_ex):
        list_lo = list_lo_theta
    else:
        list_lo = list_lo_ex


    return list_hi, list_lo



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def run_GProfiler(genelist, fileOUT):
    """
    run GProfiler GO analysis webserver through their API
    fileOUT: filename for saving the results dataframe
    """

    gp = GProfiler(return_dataframe=True)
    result = gp.profile(organism='scerevisiae', query=genelist)
    result.to_csv(fileOUT, header=True, index=False, sep='\t')

    return result




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def compile_significants(redundant_df, DT_df, DTmod_hi_df, DTmod_lo_df, DTnmod_hi_df, DTnmod_lo_df):
    """
    set a lower pvalue a threshold to focus on strong hits 
    merge all GO terms that have at least one hit at this level
    compile joint dataframe for easy visualization
    """

    p = 0.0001

    list_resultDFs = [redundant_df, DT_df, DTmod_hi_df, DTmod_lo_df, DTnmod_hi_df, DTnmod_lo_df]
    list_GO = []

    combinedDF = pd.DataFrame(columns=['source', 'GO', 'name', 'redundant', 'DT', 'DTmod_hi', 'DTmod_lo', 'DTnmod_hi', 'DTnmod_lo'])

    for GOdataframe in list_resultDFs:
        for entry in range(len(GOdataframe)):
            current_source = GOdataframe.loc[entry]['source']
            current_pval = GOdataframe.loc[entry]['p_value']
            if "GO" in current_source and current_pval < p:
                list_GO.append(GOdataframe.loc[entry]['native'])

    list_GO = list(set(list_GO))

    for GO in list_GO:
        current_pval = []
        current_source = None
        current_name = None 
        for GOdataframe in list_resultDFs:
            current_entry = GOdataframe[GOdataframe['native']==GO]
            if len(current_entry) > 0:
                current_pval.append( float(current_entry['p_value']) )
                current_source = current_entry['source'].item()
                current_name = current_entry['name'].item()
                current_name = str(current_name)
                #manually shorten some GO term names for plotting
                if current_name == "maturation of SSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)":
                    current_name = "maturation of SSU-rRNA from tricistronic rRNA transcript"
                elif current_name == "nuclear outer membrane-endoplasmic reticulum membrane network":
                    current_name = "nuclear outer membrane - ER membrane network"
            else:
                current_pval.append( np.nan)

        current_row = list([current_source, GO, current_name])
        for pval in current_pval:
            current_row.append(pval)
        combinedDF.loc[len(combinedDF)] = (current_row)

    combinedDF.sort_values('redundant', ascending=True, inplace=True)
    combinedDF.sort_values('DT', ascending=True, inplace=True)
    combinedDF = combinedDF[::-1]

    combinedDF.to_csv("../data/figures/figure6/GO_combined.txt", header=True, index=False, sep='\t', na_rep='NA')





if __name__ == '__main__':

    # load data
    mm_consensus = pickle.load(open("../data/processed/mm_consensus.pkl", 'rb'))
    sequence_codons = pickle.load(open('../data/processed/yeast_codons.pkl', 'rb'))

    # load DT sequence data -----------------------------------------
    kmers_all  = pd.read_csv("../data/figures/figure3/kmer_all.txt", header=0, index_col=False, sep='\t')
    kmers_DT   = pd.read_csv("../data/figures/figure3/kmer_filtered.txt", header=0, index_col=False, sep='\t')      
    DT_w, DT_wo = split_kmerdf(kmers_DT, ['GAA']) #['AAA', 'CAA', 'GAA']

    # get counts
    counts_all = kmer_density(kmers_all, sequence_codons, mm_consensus)
    counts_DT = kmer_density(kmers_DT, sequence_codons, mm_consensus)
    counts_DT_w = kmer_density(DT_w, sequence_codons, mm_consensus)
    counts_DT_wo = kmer_density(DT_wo, sequence_codons, mm_consensus)
  
    # compile gene lists based on tails of enrichment (density or difference)
    list_redundant_enriched, list_redundant_depleted  = get_genelist(counts_all, 'redundant', 'length', method='division')
    list_DT_enriched, list_DT_depleted = get_genelist(counts_DT, 'total', 'length', method='division')
    list_DTmod_high, list_DTmod_low = get_genelist(counts_DT_w, 'high', 'low')
    list_DTnmod_high, list_DTnmod_low = get_genelist(counts_DT_wo, 'high', 'low')


    #run GO analysis through webserver API and combine results
    result_redundant_enriched  = run_GProfiler(list_redundant_enriched, "../data/figures/figure6/GO_redundant_enriched.txt")
    result_DT_enriched         = run_GProfiler(list_DT_enriched, "../data/figures/figure6/GO_DT_enriched.txt")
    result_DTmod_high          = run_GProfiler(list_DTmod_high, "../data/figures/figure6/GO_DTmod_high.txt")
    result_DTmod_low           = run_GProfiler(list_DTmod_low, "../data/figures/figure6/GO_DTmod_low.txt")
    result_DTnmod_high         = run_GProfiler(list_DTnmod_high, "../data/figures/figure6/GO_DTnmod_high.txt")
    result_DTnmod_low          = run_GProfiler(list_DTnmod_low, "../data/figures/figure6/GO_DTnmod_low.txt")


    compile_significants(result_redundant_enriched, result_DT_enriched, result_DTmod_high, result_DTmod_low, result_DTnmod_high, result_DTnmod_low)
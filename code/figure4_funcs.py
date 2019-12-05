import os, sys
import numpy as np
import pandas as pd
import subprocess
import glob
import csv
import pickle
from scipy import stats
from Bio.Seq import Seq
from Bio import SeqIO




#--------------------------------------------------------------------------------------------
def parse_domain_boundaries(fileOUT):
    """
    parse dictionary of domain boundaries from threadom domain predictions
    alt switch: malstrom/baker superfamily paper!
    key: ORF name
    value: vector of domain boundaries
    """

    db_dict = {}

    domainTableDump = pd.read_csv("../data/accessory/domainTableDump.txt", header=0, index_col=False, sep='\t') 
    domainRegionTableDump = pd.read_csv("../data/accessory/domainRegionTableDump.txt", header=0, index_col=False, sep='\t')

    matchid = {}
    for line in open("../data/accessory/match.id", 'r'):
        current_line = line.split()
        current_id = str(current_line[0])
        current_orf = current_line[1]
        matchid[current_id] = current_orf

    list_genes = list(set(domainTableDump['geneID']))

    for gene in list_genes:
        current_orf = matchid.get(str(gene), 'none')
        current_domains = domainTableDump[domainTableDump['geneID']==gene]['domainID'].values
        current_domainboundaries = []

        for domain in current_domains:
            current_region = domainRegionTableDump[domainRegionTableDump['domainID']==domain]
            if len(current_region) == 1:
                current_start = current_region['start'].item()
                current_end = current_region['end'].item()
            elif len(current_region) > 1:
                current_start = current_region['start'].values[0]
                current_end = current_region['end'].values[-1]

            current_domainboundaries.append(current_end)

        if len(current_domainboundaries) > 1 and current_orf != 'none':     
            db_dict[current_orf] = sorted( current_domainboundaries[:-1] )
        else:
            db_dict[current_orf] = "none"

    pickle.dump(db_dict, open(fileOUT, 'wb'))

    return db_dict



#--------------------------------------------------------------------------------------------
def parse_RNAfolding(fileIN, fileOUT):
    """
    load Segal paper data into dictionary
    he's got weird reference sequences that don't match yeast codon seqs, so this function includes a
    shifty alignment (no-gaps alignment) to identify the parts that need to be chopped
    """

    parsRNA_seq = SeqIO.to_dict(SeqIO.parse('../data/accessory/RNA/sce_genes.fasta', "fasta"))
    reference_seq = SeqIO.to_dict(SeqIO.parse('../data/reference/orf_coding.fasta', "fasta"))

    rnafold_dict = {}
    rnafold_coords = {}

    weird_stuff = 0

    for line in open(fileIN, 'r'):
        current_line = line.split()
        current_orf = current_line[0]
        current_length = current_line[1]
        current_profile = np.array( current_line[2].split(';') )

        if current_orf in list(parsRNA_seq.keys()) and current_orf in list(reference_seq.keys()):
            seq_pars = np.array( list( parsRNA_seq[current_orf] ) )
            seq_ref  = np.array( list( reference_seq[current_orf] ) )
      
            if len(seq_pars) > len(seq_ref):                
                slack = len(seq_pars) - len(seq_ref) 
                aln = np.zeros(( slack ))

                for shift in range( slack ):
                    seq_test = seq_pars[shift:(shift+len(seq_ref))]
                    score = np.sum( seq_test == seq_ref )
                    aln[shift] = score

                coords_start = np.argmax(aln)    
                coords_end = coords_start + len(seq_ref) 
                coords = [coords_start, coords_end]
                current_profile_trimmed = np.array(current_profile[coords_start:coords_end], dtype=float)

                # convert from nt to codon/aa: average over 3 nts per codon/aa
                current_L = int(np.floor( len(current_profile_trimmed)/3))
                current_rnafold = np.reshape(current_profile_trimmed[0:3*current_L], ( current_L, 3) )
                current_rnafold = np.sum( current_rnafold, 1)/3.

                rnafold_coords[current_orf] = coords
                rnafold_dict[current_orf] = current_rnafold
        else:
            #print(current_orf, "what is going on here? ah, you're non-coding ... ")
            weird_stuff += 1

    pickle.dump(rnafold_dict, open(fileOUT, 'wb') ) 
 
    return rnafold_dict





#--------------------------------------------------------------------------------------------
def parse_charges():
    """
    compile dictionary of clusters of positive net charge in AA seqs
    cluster definition following C Charneski & L Hurst:
    - two positively charged residues within five amino acids, 
    - three positively charged residues within eight amino acids, 
    - four or five positively charged amino acids within 10 amino acids, 
    - and six or more positive charges within 16 amino acids,  
    """

    aa_seq = pickle.load(open("../data/processed/yeast_aa.pkl", "rb"))

    pos_charges = ['H', 'K', 'R']    
    neg_charges = ['D', 'E']

    result = {}

    cluster_dict = {5:2, 8:3, 10:5, 16:6}
    cluster_size = 8
    cluster_N = cluster_dict[cluster_size]
 
    list_orfs = list( aa_seq.keys() )

    for orf in list_orfs:
        current_seq = np.array(aa_seq[orf])
        current_clusters = []
        for pos in range( len(current_seq) - cluster_size):
            current_window = current_seq[pos:pos+cluster_size]
            score_pos = 0
            score_neg = 0
            for aa in current_window:
                if aa in pos_charges:
                    score_pos += 1
                elif aa in neg_charges:
                    score_neg += 1
            if score_pos >= cluster_N and score_neg <= 0:
                current_clusters.append(pos) 
        
        result[orf] = current_clusters

    return result





#--------------------------------------------------------------------------------------------
def parse_noptclusters():
    """
    get dictionary with positions of clusters of nonoptimal codons
    nonoptimal (following Charneski&Hurst based on bottom% of tAI): [CGA, ATA, CTT, CTG, CTC, CGG, AGT, CCC, GCG, AGC, CCT, TCG, TGT, ACG, and GTG]
    this is also the same used in Fig2 (bottom 25% of tAI)
    """

    codon_seq = pickle.load(open("../data/processed/yeast_codons.pkl", "rb"))

    # standard bottom 25% of tAI scale
    nonoptimal = ['CGA', 'ATA', 'CTT', 'CTG', 'CTC', 'CGG', 'AGT', 'CCC', 'GCG', 'AGC', 'CCT', 'TCG', 'TGT', 'ACG', 'GTG']

    result = {}

    cluster_dict = {5:2, 8:3, 10:5, 16:6}
    cluster_size = 5
    cluster_N = cluster_dict[cluster_size]

    list_orfs = list( codon_seq.keys() )

    for orf in list_orfs:
        current_seq = np.array(codon_seq[orf])
        current_clusters = []        
        for pos in range( len(current_seq) - cluster_size):
            current_window = list( current_seq[pos:pos+cluster_size] )
            score = 0
            for nopt in nonoptimal:
                if nopt in current_window:
                    score += 1
            if score >= cluster_N:
                current_clusters.append(pos) 
 
        result[orf] = current_clusters

    return result





#--------------------------------------------------------------------------------------------
def parse_loqate(fileIN='../data/accessory/loqate.txt'):
    """
    parse localisation data
    """  
    loqate_dict = {}

    for line in open(fileIN, 'r'):
        current_line = line.rstrip('\n').split()
        current_orf = current_line[0]
        current_loc = current_line[1]
        loqate_dict[current_orf] = current_loc

    return loqate_dict



#--------------------------------------------------------------------------------------------
def test_kmer(kmerdf):
    """
    compile DF with association stats for each kmer in input data
    class is 1 for above threshold RD and 0 for below threshold RD
    feature is 1 for present and 0 for absent
    -> OR > 1 means an enrichment of high RD associated with feature
    -> OR < 1 means depletion of high RD associated with feawture
    """
 
    theta_tunnel = 35 				

    result_df = pd.DataFrame(columns=['feature', 'OR', 'pval'])

    table_rnafold = np.zeros(( 2, 2 ))
    table_charges = np.zeros(( 2, 2 ))
    table_nopt    = np.zeros(( 2, 2 ))
    table_loc_ER   = np.zeros(( 2, 2 ))
    table_loc_mito = np.zeros(( 2, 2 ))
    table_loc_nuc  = np.zeros(( 2, 2 ))

    for i in range( len(kmerdf) ):
        current_kmer = kmerdf.iloc[i]['kmer']
        current_orf  = kmerdf.iloc[i]['ORF']
        current_position = kmerdf.iloc[i]['position']
        current_class = int( kmerdf.iloc[i]['class'] )


        # rna folding strength --------------------------------------------
        current_rnafold = rnafold.get(current_orf, np.array([np.nan]))   
        if len(current_rnafold) > 1 and current_position > 50:				# if no RNA fold data availalbe for gene, it's set to [np.nan]
            current_rnafold_atoffset = current_rnafold[(current_position+4)]	#  -12nt, so -4codons, -3codons may be too close, but downstream should have more effect than upstream
            #if np.any(current_rnafold_atoffset) > 1:
            if current_rnafold_atoffset > 0:
                current_feature_rnafold = 1
            else:
                current_feature_rnafold = 0        
            table_rnafold[current_feature_rnafold, current_class] += 1

        # clusters of positive charges -------------------------------------
        current_chargeclust = chargedclusters[current_orf]
        if len(current_chargeclust) > 0 :
            current_dist_charges = current_position - current_chargeclust 
            current_dist_charges = current_dist_charges[current_dist_charges >= 0 ] 		# only inside the tunnel but after constriction site
            if np.any( current_dist_charges < theta_tunnel ):
                current_feature_charges = 1
            else:
                current_feature_charges = 0      
            table_charges[current_feature_charges, current_class] += 1

        # clusters of nonoptimal codons ------------------------------------
        current_noptclust = noptclusters[current_orf]
        if len(current_noptclust) > 0 :
            current_dist_nopt = current_position - current_noptclust
            current_dist_nopt = current_dist_nopt[current_dist_nopt <= 0] 	# only what was translated just before
            if np.any( current_dist_nopt > - 3 ):							# just before, maybe even smaller is more meaningful
                current_feature_nopt = 1
            else:
                current_feature_nopt = 0
            table_nopt[current_feature_nopt, current_class] += 1

        # loqate localisation-----------------------------------------------
        current_loc = loqate.get(current_orf, 'none' )
        #if current_loc != 'below_threshold' and current_loc != 'none':
        if current_loc == 'cytosol' or current_loc == 'ER':
            if 'cyto' in current_loc: # == 'cytosol':
                current_feature_ER = 0
            else:
                current_feature_ER = 1
            table_loc_ER[current_feature_ER, current_class] += 1

        if current_loc == 'cytosol' or current_loc == 'mitochondria':
            if 'cyto' in current_loc: # == 'cytosol':
                current_feature_mito = 0
            else:
                current_feature_mito = 1
            table_loc_mito[current_feature_mito, current_class] += 1

        if current_loc == 'cytosol' or current_loc == 'nucleus':
            if 'cyto' in current_loc: # == 'cytosol':
                current_feature_nuc = 0
            else:
                current_feature_nuc = 1
            table_loc_nuc[current_feature_nuc, current_class] += 1


  
    OR_rnafold, pval_rnafold = stats.fisher_exact(table_rnafold)
    OR_charges, pval_charges = stats.fisher_exact(table_charges)
    OR_nopt, pval_nopt = stats.fisher_exact(table_nopt)
   
    OR_loc_ER, pval_loc_ER = stats.fisher_exact(table_loc_ER)
    OR_loc_mito, pval_loc_mito = stats.fisher_exact(table_loc_mito)
    OR_loc_nuc, pval_loc_nuc = stats.fisher_exact(table_loc_nuc)


    result_df.loc[len(result_df)] = ('nopt', OR_nopt, pval_nopt)
    result_df.loc[len(result_df)] = ('rnafold', OR_rnafold, pval_rnafold)
    result_df.loc[len(result_df)] = ('poscharges', OR_charges, pval_charges)
    result_df.loc[len(result_df)] = ('locER', OR_loc_ER, pval_loc_ER)
    result_df.loc[len(result_df)] = ('locMito', OR_loc_mito, pval_loc_mito)
    result_df.loc[len(result_df)] = ('locNuke', OR_loc_nuc, pval_loc_nuc)

   
    return result_df



#--------------------------------------------------------------------------------------------
def test_positionbias(kmertable, data_mc, data_mm):
    """
    check if there is any positional bias
    RL: relative length; not good but else difficult to compare seqs of diff length
    """

    reference_seq = SeqIO.to_dict(SeqIO.parse('../data/reference/orf_coding.fasta', "fasta"))

    resultDF = pd.DataFrame(columns=['class', 'position', 'mc'])
    for i in range(len(kmertable)):
        current_orf = kmertable.iloc[i]['ORF']
        current_pos = kmertable.iloc[i]['position']
        current_class = kmertable.iloc[i]['class']
        current_mc = np.mean( data_mc[current_orf][current_pos:(current_pos+3)] )

        current_len = np.floor( len(reference_seq[current_orf]) / 3. )  
     
        if current_class == 1:
            resultDF.loc[len(resultDF)] = ("pos"+str(current_class), current_pos, current_mc)
   
    resultDF.to_csv("../data/figures/figure4/DT_clusters.txt", header=True, index=False, sep='\t')


    def bg_positionbias(data_mc, data_mm):
        """
        check if there is a general position bias in the mc data
        """

        theta = pd.read_csv("../data/figures/figure3/theta.txt", header=0, index_col=False, sep='\t')

        trim = 20			# omit first and last 20 codons per gene due to known biases of translation initiation and termination
        kmer = 3

        resultDF = pd.DataFrame(columns=['class', 'position', 'mc'])

        list_orfs = list( data_mc.keys() )

        for ix, orf in enumerate( list_orfs ) :
            print(ix, orf)
            current_consensus = data_mc[orf]
            current_mm = data_mm[orf]

            current_theta_lo10 =  theta[theta['ORF']==orf]['p3_10'].item()
            current_theta_hi90 =  theta[theta['ORF']==orf]['p3_90'].item()    

            for pos in range( trim, len(current_consensus) - (trim + kmer) ):	# omit first/last 20 positions and allow for kmer length
                current_score = np.mean(current_consensus[pos:pos+kmer])
                current_pos = pos #/ len(current_consensus)
 
                if current_score > current_theta_hi90:
                    resultDF.loc[len(resultDF)] = ("bg1", current_pos, current_score) 

                #elif current_score <= current_theta_lo10:
                #    resultDF.loc[len(resultDF)] = ("bg0", current_pos, current_score)

        
        resultDF.to_csv("../data/figures/figure4/bg_clusters.txt", header=True, index=False, sep='\t')


    bg_positionbias(data_mc, data_mm)




#--------------------------------------------------------------------------------------------
def test_clusterbias(kmertable, data_mc, data_mm):
    """
    test if there is a bias in broader peaks (bleed-over of signal) vs. sharp peaks
    """

    def kmer_clusterbias(kmertable, data_mc, data_mm):
        """
        subfunction to check if there is a general position bias in the mc data
        """

        theta = pd.read_csv("../data/figures/figure3/theta.txt", header=0, index_col=False, sep='\t')

        trim = 20			# omit first and last 20 codons per gene due to known biases of translation initiation and termination
        kmer = 3
        w = 3				# 3 to either side, total window size of 7

        resultDF = pd.DataFrame(columns=['class', 'cluster', 'position', 'mc'])

        list_orfs = list( data_mc.keys() )

        tripletdict = {}

        n_class1_clust = 0
        n_class1_nclust = 0

        for ix, orf in enumerate( list_orfs ) :
            current_consensus = data_mc[orf]
            current_mm = data_mm[orf]
            current_consensus[current_mm == False] = np.nan 

            current_theta_lo10 =  theta[theta['ORF']==orf]['p3_10'].item()
            current_theta_hi90 =  theta[theta['ORF']==orf]['p3_90'].item()    

            score = np.zeros(( len(current_consensus) )) * np.nan 

            for pos in range( trim, len(current_consensus) - (trim + kmer) ):	# omit first/last 20 positions and allow for kmer length
                score[pos] = np.mean(current_consensus[pos:pos+kmer])

            tripletdict[orf] = score

        for i in range(len(kmertable)):
            current_orf = kmertable.iloc[i]['ORF']
            current_pos = kmertable.iloc[i]['position']
            current_class = kmertable.iloc[i]['class']
            current_mc = np.mean( data_mc[current_orf][current_pos:(current_pos+3)] )

            score = tripletdict[current_orf]

            if current_mc > current_theta_hi90:
                current_window = score[(current_pos-w):(current_pos+w)]
                current_window = current_window[~np.isnan(current_window)]
                if np.sum( current_window > current_theta_hi90 ) >= 3:
                    resultDF.loc[len(resultDF)] = ("DT_cl", 1, current_pos, current_mc)
                    n_class1_clust += 1
                else:
                    resultDF.loc[len(resultDF)] = ("DT_n", 0, current_pos, current_mc)
                    n_class1_nclust += 1


        resultDF.to_csv("../data/figures/figure4/cluster_kmer_position.txt", header=True, index=False, sep='\t')
    
        return n_class1_clust, n_class1_nclust, resultDF



    def bg_clusterbias(data_mc, data_mm):
        """
        subfunction to check if there is a general position bias in the mc data
        for speed only look at high RD clusters, as the low ones may include missing coverage
        """

        theta = pd.read_csv("../data/figures/figure3/theta.txt", header=0, index_col=False, sep='\t')

        trim = 20			# omit first and last 20 codons per gene due to known biases of translation initiation and termination
        kmer = 3
        w = 3				# 3 to either side, total window size of 7

        resultDF = pd.DataFrame(columns=['class', 'cluster', 'position', 'mc'])

        list_orfs = list( data_mc.keys() )

        n_class1_clust = 0
        n_class1_nclust = 0

        for ix, orf in enumerate( list_orfs ) :
            print(ix, orf)
            current_consensus = data_mc[orf]
            current_mm = data_mm[orf]
            current_consensus[current_mm == False] = np.nan 

            current_theta_lo10 =  theta[theta['ORF']==orf]['p3_10'].item()
            current_theta_hi90 =  theta[theta['ORF']==orf]['p3_90'].item()    

            score = np.zeros(( len(current_consensus) )) * np.nan 

            for pos in range( trim, len(current_consensus) - (trim + kmer) ):	# omit first/last 20 positions and allow for kmer length
                score[pos] = np.mean(current_consensus[pos:pos+kmer])
                current_pos = pos #/ len(current_consensus)
 
            for pos in range( trim, len(current_consensus) - (trim + kmer) ):

                if score[pos] > current_theta_hi90:
                    current_window = score[(current_pos-w):(current_pos+w)]
                    current_window = current_window[~np.isnan(current_window)]
                    if np.sum( current_window > current_theta_hi90 ) >= 3:
                        resultDF.loc[len(resultDF)] = ("BG_cl", 1, pos, score[pos])
                        n_class1_clust += 1
                    else:
                        resultDF.loc[len(resultDF)] = ("BG_n", 0, pos, score[pos])
                        n_class1_nclust += 1

        resultDF.to_csv("../data/figures/figure4/cluster_bg_position.txt", header=True, index=False, sep='\t')

        return n_class1_clust, n_class1_nclust, resultDF



    kmer_n_clust, kmer_n_nclust, kmer_DF = kmer_clusterbias(kmertable, data_mc, data_mm)
    bg_n_clust, bg_n_nclust, bg_DF = bg_clusterbias(data_mc, data_mm)


    resultDF = pd.DataFrame(columns=['category', 'cluster', 'not'])
    resultDF.loc[len(resultDF)] = ("DT", kmer_n_clust, kmer_n_nclust )
    resultDF.loc[len(resultDF)] = ("BG", bg_n_clust, bg_n_nclust )

    resultDF.to_csv("../data/figures/figure4/clusters.txt", header=True, index=False, sep='\t')






if __name__ == '__main__':


    # load auxiliary data -------------------------------------------

    # if pickled dictionary exists load it, else generate it
    PATH_domainboundaries = "../data/accessory/domainboundaries.pkl"
    if os.path.exists(PATH_domainboundaries):
        domainboundaries = pickle.load(open(PATH_domainboundaries, "rb"))
    else:
        domainboundaries = parse_domain_boundaries(PATH_domainboundaries)

    # if pickled dictionary exists load it, else generate it
    PATH_RNAfold = "../data/accessory/rnafold.pkl"
    if os.path.exists(PATH_RNAfold):
        rnafold = pickle.load(open(PATH_RNAfold, "rb"))
    else:
        rna_file = '../data/accessory/RNA/sce_Score.tab'
        rnafold = parse_RNAfolding(rna_file, PATH_RNAfold)

    ssb_binding = pd.read_csv("../data/accessory/ssb_biding_split_df.csv")
    chargedclusters = parse_charges()
    noptclusters = parse_noptclusters()
    loqate = parse_loqate('../data/accessory/loqate.txt')



    scikit_consensus = pickle.load(open("../data/processed/mc_dict.pkl", 'rb'))
    mm_consensus = pickle.load(open("../data/processed/mm_consensus.pkl", 'rb'))




    # load DT sequence data -----------------------------------------
    kmers = pd.read_csv("../data/figures/figure3/kmer_all.txt", header=0, index_col=False, sep='\t')
    kmers_DT   = pd.read_csv("../data/figures/figure3/kmer_filtered.txt", header=0, index_col=False, sep='\t')
    kmers_nDT_lo = pd.read_csv("../data/figures/figure3/kmer_filtered_nonDT-.txt", header=0, index_col=False, sep='\t')
    kmers_nDT_hi = pd.read_csv("../data/figures/figure3/kmer_filtered_nonDT+.txt", header=0, index_col=False, sep='\t')



    test_positionbias(kmers_DT, scikit_consensus, mm_consensus)
    test_clusterbias(kmers_DT, scikit_consensus, mm_consensus)


    print("Analysis of DT sequences")
    output_DT = test_kmer(kmers_DT)
    output_DT.to_csv("../data/figures/figure4/associations_kmers_DT.txt", header=True, index=False, sep='\t')
    print(output_DT)


    print("Analysis of nDT-")
    output_nDT_lo = test_kmer(kmers_nDT_lo)
    output_nDT_lo.to_csv("../data/figures/figure4/associations_kmers_nDT-.txt", header=True, index=False, sep='\t')
    print(output_nDT_lo)


    print("Analysis of nDT+")
    output_nDT_hi = test_kmer(kmers_nDT_hi)
    output_nDT_hi.to_csv("../data/figures/figure4/associations_kmers_nDT+.txt", header=True, index=False, sep='\t')
    print(output_nDT_hi)

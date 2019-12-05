import os, sys
import numpy as np
import pandas as pd
import pickle



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def test_ssb(kmertable, data_mc, data_mm, ssb):
    """
    get enrichment/density profiles of RD outside tunnel assoc'd with feature
    """
 
    L = 70
    N_rand = 100

    def local_density(data_mc, data_mm, ssb, first=False, random=False):
        """
        subfunction to switch between randomized or not
        """

        output = np.zeros(( len(ssb), L ))

        for i in range( len(ssb) ):
            current_orf = ssb.loc[i]['geneID']
            if current_orf in data_mc.keys():
                current_start = ssb.loc[i]['site_start']
                current_order = ssb.loc[i]['order']
                current_kmertable = kmertable[kmertable['ORF']==current_orf]

                current_mc = data_mc[current_orf]
                current_mm = data_mm[current_orf]
                current_mc[current_mm == False] = np.nan 
               
                if current_start < len(current_mc) - L: 

                    if random:
                        rand_mc = np.copy(current_mc[(current_start):(current_start+L)])
                        local_state = np.random.RandomState(seed=None)
                        domain_mc = local_state.choice(rand_mc, len(rand_mc), replace=False)
                    else:
                        domain_mc = current_mc[(current_start):(current_start+L)]

                    if len(domain_mc) == L:
                        if first:
                            if current_order == 1:
                                output[i,:] = domain_mc
                        else:
                            output[i,:] = domain_mc

        output = output[np.nansum(output,1)>0,:]
        output = np.nansum(output,0) / np.sum(np.isnan(output)==False,0)  

        return output



    def DT_enrichment(kmertable):
        """
        subfunction to check if region outside tunnel is enriched in sets of kmers
        """

        local_class = pd.DataFrame(columns=['class', 'position'])

        for i in range( len(ssb) ):
            current_orf = ssb.loc[i]['geneID']
            if current_orf in data_mc.keys():

                current_start = ssb.loc[i]['site_start']
                current_order = ssb.loc[i]['order']
                current_kmertable = kmertable[kmertable['ORF']==current_orf]

                current_mc = data_mc[current_orf]
                current_mm = data_mm[current_orf]
                current_mc[current_mm == False] = np.nan 

                if current_start < len(current_mc) - L: 
                    current_hi = current_kmertable[current_kmertable['class']==1]
                    current_lo = current_kmertable[current_kmertable['class']==0]

                    if len(current_hi) > 0:
                        current_hi_pos = np.array(current_hi['position']) - current_start
                        current_hi_pos = current_hi_pos[current_hi_pos > 0]		# only downstream, i.e. translated after
                        current_hi_pos = current_hi_pos[current_hi_pos < L]		# limit analysis to 70 positions downstream, i.e. ~35 for tunnel plus some
                        if len(current_hi_pos) > 0:
                            for pos in current_hi_pos:
                                local_class.loc[len(local_class)] = (1, pos)

                    if len(current_lo) > 0:
                        current_lo_pos = np.array(current_lo['position']) - current_start
                        current_lo_pos = current_lo_pos[current_lo_pos > 0]		# only downstream, i.e. translated after
                        current_lo_pos = current_lo_pos[current_lo_pos < L]		# limit analysis to 70 positions downstream, i.e. ~35 for tunnel plus some
                        if len(current_lo_pos) > 0:
                            for pos in current_lo_pos:
                                local_class.loc[len(local_class)] = (0, pos)

        return local_class


    out_rand = np.zeros(( N_rand, L ))
    for i in range(N_rand):
        out_rand[i,:] = local_density(data_mc, data_mm, ssb, first=True, random=True)
    rand_mean = np.mean(out_rand, 0)
    rand_sd   = np.std(out_rand, 0) 

    result_onlyfirst = local_density(data_mc, data_mm, ssb, first=True, random=False)
    result_all = local_density(data_mc, data_mm, ssb, first=False, random=False)

    resultDF = pd.DataFrame({'position': list(np.arange(len(result_onlyfirst))),
    	                    'firstssb': list(result_onlyfirst), 
                            'allssb': list(result_all), 
                            'rand_mean': list(rand_mean), 
                            'rand_sd': list(rand_sd) } )  
    resultDF.to_csv("../data/figures/figure5/ssb_enrichment.txt", header=True, index=False, sep='\t')

    
    DTkmers_byclass = DT_enrichment(kmertable)
    DTkmers_byclass.to_csv("../data/figures/figure5/ssb.DTkmers", header=True, index=False, sep='\t')



   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def test_domains(kmertable, data_mc, data_mm, domains, ssb):
    """
    get enrichment/density profiles of RD outside tunnel of domainboundaries
    first: test only first internal domain boundary or all
    chap: test only domain boundaries that are placed before any SSB binding site
    random: randomize consensus RD values in downstream window of length L
    """

    L = 70
    N_rand = 100


    def local_density(data_mc, data_mm, domains, ssb, first=False, chap=False, random=False):
        """
        subfunction to switch between randomized or not
        random: boolean whether to scramble consensus RD or not
        first: boolean whether to only consider the first domain boundary or all
        """

        output = np.zeros(( 30000, L ))			# initialize with some large value
        i = 0
        for orf in domains.keys():
            if orf in data_mc.keys():
                current_domains = domains[orf]
                current_mc = data_mc[orf]
                current_mm = data_mm[orf]
                current_mc[current_mm == False] = np.nan 
            
                if chap: 									# only domain boundaries before any SSB binding site
                    current_ssb = ssb[ssb['geneID']==orf]
                    current_ssb_pos = np.array(current_ssb['site_start'])
                    filtered_domains = []
                    if current_domains != 'none' and len(current_ssb_pos) > 0:
                        for domain in list(current_domains):
                            current_dists = int(domain) - current_ssb_pos
                            if np.all(current_dists > 0):
                                filtered_domains.append(domain)    
                        if len(filtered_domains) > 0:
                            current_domains = np.copy(np.array(filtered_domains))
                        else:
                            current_domains = "none"
                    elif current_domains != "none" and len(current_ssb_pos) == 0:
                        current_domains = current_domains
                    else:
                        current_domains = "none" 

                if current_domains != "none":				# set to 'none' for prots without internal domain boundaries
                    current_domains = list(current_domains)

                    if first:	
                        if len(current_domains) > 1:							# limit analysis to first domain boundary if there are several
                            current_domains = [current_domains[0]]

                    for current_start in current_domains: 
                        if current_start < len(current_mc) - L: 
                            if random:
                                rand_mc = np.copy(current_mc[(current_start):(current_start+L)])
                                local_state = np.random.RandomState(seed=None)
                                domain_mc = local_state.choice(rand_mc, len(rand_mc), replace=False)
                            else:
                                domain_mc = current_mc[(current_start):(current_start+L)]

                            if len(domain_mc) == L:
                                output[i,:] = domain_mc
                                i += 1

        output = output[np.nansum(output,1)>0,:]
        output = np.nansum(output,0) / np.sum(np.isnan(output)==False, 0)  # average profile ignoring nan's, e.g. those from the MM positions

        return output


    def DT_enrichment(kmertable):
        """
        subfunction to check if region outside tunnel is enriched in sets of kmers
        """

        local_class = pd.DataFrame(columns=['class', 'position'])

        for orf in domains.keys():
            if orf in data_mc.keys():
                current_domains = domains[orf]
                current_mc = data_mc[orf]
                current_mm = data_mm[orf]
                current_mc[current_mm == False] = np.nan 
                current_kmertable = kmertable[kmertable['ORF']==orf]

                if current_domains != "none":				# set to 'none' for prots without internal domain boundaries
                    current_domains = list(current_domains)           
                    for current_start in current_domains: 
                        if current_start < len(current_mc) - L: 
                            current_hi = current_kmertable[current_kmertable['class']==1]
                            current_lo = current_kmertable[current_kmertable['class']==0]

                            if len(current_hi) > 0:
                                current_hi_pos = np.array(current_hi['position']) - current_start
                                current_hi_pos = current_hi_pos[current_hi_pos > 0]		# only downstream, i.e. translated after
                                current_hi_pos = current_hi_pos[current_hi_pos < L]	# limit analysis to 70 positions downstream, i.e. ~35 for tunnel plus some
                                if len(current_hi_pos) > 0:
                                    for pos in current_hi_pos:
                                        local_class.loc[len(local_class)] = (1, pos)

                            if len(current_lo) > 0:
                                current_lo_pos = np.array(current_lo['position']) - current_start
                                current_lo_pos = current_lo_pos[current_lo_pos > 0]		# only downstream, i.e. translated after
                                current_lo_pos = current_lo_pos[current_lo_pos < L]		# limit analysis to 70 positions downstream, i.e. ~35 for tunnel plus some
                                if len(current_lo_pos) > 0:
                                    for pos in current_lo_pos:
                                        local_class.loc[len(local_class)] = (0, pos)

        return local_class


    result_first = local_density(data_mc, data_mm, domains, ssb, first=True, chap=False, random=False)
    result_all = local_density(data_mc, data_mm, domains, ssb, first=False, chap=False, random=False)
    result_nonssb = local_density(data_mc, data_mm, domains, ssb, first=True, chap=True, random=False)
    
    out_rand = np.zeros(( N_rand, L ))
    for i in range(N_rand):
        out_rand[i,:] = local_density(data_mc, data_mm, domains, ssb, first=True, chap=False, random=True)
    rand_mean = np.mean(out_rand, 0)
    rand_sd   = np.std(out_rand, 0) 

    resultDF = pd.DataFrame({'position': list(np.arange(len(result_first))),
    	                    'first_db': list(result_first), 
    	                    'all_db': list(result_all),
    	                    'nonssb': list(result_nonssb),
                            'rand_mean': list(rand_mean), 
                            'rand_sd': list(rand_sd) } )  
    resultDF.to_csv("../data/figures/figure5/domain_RD.txt", header=True, index=False, sep='\t')

    
    DTkmers_byclass = DT_enrichment(kmertable)
    DTkmers_byclass.to_csv("../data/figures/figure5/domain.DTkmers", header=True, index=False, sep='\t')






if __name__ == '__main__':


    # load base data ------------------------------------------------
    scikit_consensus = pickle.load(open("../data/processed/mc_dict.pkl", 'rb'))
    mm_consensus = pickle.load(open("../data/processed/mm_consensus.pkl", 'rb'))

    # load auxiliary data -------------------------------------------
    domainboundaries = pickle.load(open("../data/accessory/domainboundaries.pkl", "rb"))
    ssb_binding = pd.read_csv("../data/accessory/ssb_biding_split_df.csv")
   
    # load DT sequence data -----------------------------------------
    kmers_DT   = pd.read_csv("../data/figures/figure3/kmer_filtered.txt", header=0, index_col=False, sep='\t')


    # FUNCTIONS START HERE ------------------------------------------
    test_ssb(kmers_DT, scikit_consensus, mm_consensus, ssb_binding)
    test_domains(kmers_DT, scikit_consensus, mm_consensus, domainboundaries, ssb_binding)




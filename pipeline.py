import os, sys
import numpy as np
import pandas as pd
import subprocess
import glob
import json
import csv
import pickle
import pysam
from Bio.Seq import Seq




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def download_SRA(SRA):
    """
    download sequencing file from SRA archive
    requires local install of SRA tools in path
    requires verification of filenames and paths
    """

    print("Downloading SRA archive")
    output = subprocess.run(['prefetch', '-f', 'yes', SRA], stderr=subprocess.STDOUT)

    print("Extracting FASTQ data")
    output = subprocess.run(['fastq-dump', '--gzip', NCBI_DIR+SRA+'.sra'], stderr=subprocess.STDOUT)

 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def trim_adapters(SRA):
    """
    trim sequencing adapters
    """

    cmd = 'cutadapt' + ' ' + '--cores=8' + ' ' + '--discard-untrimmed' + ' ' + '-u' + ' ' + str(1) + ' ' + '-m' + ' ' + str(10) + ' ' + '-q' + ' ' +  str(13) + ' ' + '-a' + ' ' + 'file:data/reg_adapters.fa' + ' ' + '-o' + ' ' +  TMP_DIR+SRA+'_trimmed.fastq' + ' ' + SRA+'.fastq.gz'

    print("Trimming sequencing adapters")
    output = subprocess.check_output(cmd, shell=True) #, stderr=subprocess.STDOUT).decode('UTF-8')
    


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def build_indices(genome_fasta, genome_gtf, rRNA_fasta, transcriptome_fasta):
    """
    put path to indices / pass paths as arg, e.g. STAR_DIR
    """
    
    if not os.path.exists("data/indices"):
        os.mkdir("data/indices")

  
    # 1. Bowtie index
    print("Building Bowtie index")
    if not os.path.exists(BOWTIE_DIR):
        os.mkdir(BOWTIE_DIR)
    cmd_bowtie = 'bowtie-build' + ' ' + genome_fasta + ' ' + BOWTIE_DIR+'/yeast'
    output = subprocess.run(cmd_bowtie, shell=True)

    cmd_rRNA = 'bowtie-build' + ' ' + rRNA_fasta + ' ' + BOWTIE_DIR+'/rRNA'
    output = subprocess.run(cmd_rRNA, shell=True)
    
    # 2. STAR index
    print("Building STAR index")
    if not os.path.exists(STAR_DIR):
        os.mkdir(STAR_DIR)
    cmd_STAR = 'STAR' + ' ' + '--runThreadN' + ' ' + '4' + ' ' + '--runMode' + ' ' + 'genomeGenerate' + ' ' + '--genomeDir' + ' ' + STAR_DIR + ' ' + '--genomeFastaFiles' + ' ' + genome_fasta + ' ' + '--sjdbGTFfile' + ' ' + genome_gtf #+ ' ' + '--sjdbOverhang' + ' ' + 'max(ReadLength)-1'
    output = subprocess.run(cmd_STAR, shell=True)


#    run build transcriptome fasta. 
    if not os.path.exists(STAR_TRANSCRIPTOME_DIR):
        os.mkdir(STAR_TRANSCRIPTOME_DIR)
    cmd_STAR = 'STAR' + ' ' + '--runThreadN' + ' ' + '4' + ' ' + '--runMode' + ' ' + 'genomeGenerate' + ' ' + '--genomeDir' + ' ' + STAR_TRANSCRIPTOME_DIR + ' ' + '--genomeFastaFiles' + ' ' + transcriptome_fasta # + ' ' + '--sjdbGTFfile' + ' ' + genome_gtf #+ ' ' + '--sjdbOverhang' + ' ' + 'max(ReadLength)-1'
    output = subprocess.run(cmd_STAR, shell=True)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def map_reads(SRA):
    """
    maps reads (bowtie to rRNA for legacy?) to extract ambiguous and uniquely mapped reads
    """

    #1. bowtie to rRNA
    print("Bowtie alignement on contaminant RNA...")
    cmd_bowtie = 'bowtie'+ ' ' + '-a' + ' ' + '-p6' + ' ' + '-S' + ' ' + '--un' + ' ' +  TMP_DIR+SRA+'_rrnaUnmapped.fastq' + ' ' + BOWTIE_DIR+'/rRNA' + ' ' + TMP_DIR+SRA+'_trimmed.fastq' + ' ' + '|' + ' ' + 'samtools view -@ 6 -bS' + ' ' + '>' + TMP_DIR+SRA+'_trimmed_rrna.bam'
    output = subprocess.run(cmd_bowtie, shell=True)

    # 2. STAR to ref genome
    print("STAR alignement to yeast genome...")
    cmd_STAR = 'STAR --outSAMtype BAM Unsorted --runThreadN 6 --winAnchorMultimapNmax 200 --seedSearchStartLmax 15 --genomeDir' + ' ' + STAR_DIR + ' ' + '--readFilesIn' + ' ' + TMP_DIR+SRA+'_rrnaUnmapped.fastq' + ' ' + '--outFileNamePrefix' + ' ' + TMP_DIR+SRA+'_STAR_'
    output = subprocess.run(cmd_STAR, shell=True)

    # 3. Samtools keep uniquely mapped reads and sort
    print("Samtools to keep uniquely mapped reads and sort...")
    cmd_samtools1 = 'samtools view -@ 6 -b -q 255 -o' + ' ' + TMP_DIR+SRA+'_yeast_uniqueMapped_reads.bam' + ' ' + TMP_DIR+SRA+'_STAR_Aligned.out.bam'
    output = subprocess.run(cmd_samtools1, shell=True)

    cmd_samtools2 = 'samtools sort -@ 6 -o' + ' ' + TMP_DIR+SRA+'_yeast_uniqueMapped_reads_sorted.bam' + ' ' + TMP_DIR+SRA+'_yeast_uniqueMapped_reads.bam'
    output = subprocess.run(cmd_samtools2, shell=True)

    cmd_samtools3 = 'samtools index' + ' ' + TMP_DIR+SRA+'_yeast_uniqueMapped_reads_sorted.bam'
    output = subprocess.run(cmd_samtools3, shell=True)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def run_scikit_ribo(SRA, genome_fasta, genome_gtf):
    """
    wrapper to run scikit-ribo from the same pipeline
    requires local install of modified scikit-ribo toolbox
    requires local install of all dependencies of scikit-ribo environment (see conda environment file)
    """

    # 3. Scikit-ribo index
    print("Building scikit-ribo index")
    if not os.path.exists(SCIKIT_DIR):
        os.mkdir(SCIKIT_DIR)
    cmd_scikit = 'python' + ' ' + SCIKIT_PATH + 'scikit-ribo-build.py' + ' ' + '-g' + ' ' + genome_gtf + ' ' + '-f' + ' ' + genome_fasta + ' ' + '-p' + ' ' + SRA + ' ' + '-o' + SCIKIT_DIR
    output = subprocess.run(cmd_scikit, shell=True)


    print("scikit-ribo-run.py...")
    cmd_scikit = 'python' + ' ' + SCIKIT_PATH + 'scikit-ribo-run.py' + ' ' + '-i' + ' ' + TMP_DIR+SRA+'_yeast_uniqueMapped_reads_sorted.bam' + ' ' + '-f' + ' ' + SCIKIT_DIR + ' ' + '-p' + ' ' + SRA + ' ' + '-o' + ' ' + 'TMP/scikit_'+SRA 
    output = subprocess.run(cmd_scikit, shell=True)


    print("plot_ribo_density_dict.py...")
    cmd_scikit = 'python' + ' ' + SCIKIT_PATH + 'plot_ribo_density_dict_noCDT.py' + ' ' + '-i' + ' ' + TMP_DIR+'scikit_'+SRA+'/riboseq_input.txt' + ' ' + '-g' + ' ' + 'all' + ' ' + '-o' + ' ' + TMP_DIR+'scikit_'+SRA #+'_profiles'
    output = subprocess.run(cmd_scikit, shell=True)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def run_multimapping(SRA):
    """
    identify all reads that map ambigously and their positions
    """

    if not os.path.exists("TMP/ambiguous_reads/"):
        os.mkdir("TMP/ambiguous_reads/")

    cmd_STAR = 'STAR --outSAMtype BAM SortedByCoordinate --runThreadN 8 --winAnchorMultimapNmax 200 --seedSearchStartLmax 15 --genomeDir' + ' ' + STAR_TRANSCRIPTOME_DIR + ' ' + '--readFilesIn' + ' ' + TMP_DIR+SRA+'_rrnaUnmapped.fastq' + ' ' + '--outFileNamePrefix'  + ' ' + TMP_DIR+'ambiguous_reads/'+SRA+'_STAR_transcriptome_'
    output = subprocess.run(cmd_STAR, shell=True)

 
    # Keep only multi-mapping reads:
    cmd_filter = 'python code/sam_STAR_mapq_filtering.py' + ' ' + TMP_DIR+'ambiguous_reads/'+SRA+'_STAR_transcriptome_Aligned.sortedByCoord.out.bam' + ' ' + TMP_DIR+'ambiguous_reads/'+SRA+'_STAR_transcriptome_multi_mapped_sorted.bam' + ' ' + 'all'
    output = subprocess.run(cmd_filter, shell=True)

    cmd_samtools2 = 'samtools index' + ' ' + TMP_DIR+'ambiguous_reads/'+SRA+'_STAR_transcriptome_multi_mapped_sorted.bam'
    output = subprocess.run(cmd_samtools2, shell=True)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def parse_dataframes(genome_gtf, sralist):
    """
    consolidate results into dataframes
    """

    def gather_strand_by_geneID_dict(genome_gtf):
        """
        Returns dictionary with strand orientation as values and geneIDs as Keys/
        e.g.: {'YAL012W': '+',
               'YAL069W': '+',
               'YAL068W-A': '+',
        """
        strand_by_geneID_dict = {}
        with open(genome_gtf) as f: 
            for line in f: 
                current_line = line.split('\t')
                if current_line[2] == "CDS":
                    current_orf = current_line[8].split(';')[2].split()[1].strip('\"')
                    current_strand = current_line[6]
                    strand_by_geneID_dict[current_orf] = current_strand
        return strand_by_geneID_dict


    def import_scikit_data(sralist):
        """
        Import results from scikit pipeline for all datasets contained in datsets_names.
        """
        scikit_data_dict = {}
        for dataset in sralist:
            with open(TMP_DIR+'scikit_'+dataset+'/ALL_genes_profile_dict.json', 'r') as scikit_data:
                scikit_data_dict[dataset] = [json.load(scikit_data)]
        return scikit_data_dict


    def build_mat_scikit_strandOriented(sralist, scikit_data):
        """
        Building of scikit_df based on the output of plot_ribo_density_dict.py script.

        C/-/reverse/complementary strand are taken into account and the profile values
        ("codon_density_profile", "codon_triplet", "codon_AA") are reversed. This is
        performed by adding [::-1] to C strands profile ends.

        Same profile values are also have their extremities trimmed out of 8 codons.
        (This is because the scikit-ribo pipeline considers 8 extra codons on each end,
        but here we are only interested in the coding sequence). This is performed by
        adding [8:-8] to profile lists ends.
        """

        scikit_mat = {}
        seq_codons = {}
        seq_aa = {}

        for geneID in scikit_data[sralist[0]][0].keys():
            for ix, dataset in enumerate(sralist):

                if geneID in scikit_data[dataset][0].keys():
                    current_profile = scikit_data[dataset][0].get(geneID, np.nan)
                    current_ribo = current_profile[0]
                    current_ribo = current_ribo[8:-8]
                    N = len(sralist)
                    M = len(current_ribo)
                    print(geneID, M)

                    if ix == 0:
                        current_matrix = np.zeros((N,M)) * np.nan

                        current_seq_codons = current_profile[1]
                        current_seq_codons = current_seq_codons[8:-8]

                        current_seq_aa = current_profile[2]
                        current_seq_aa = current_seq_aa[8:-8]

                        if strand_by_geneID_dict.get(geneID, "NA") == "+":
                            seq_codons[geneID] = current_seq_codons
                            seq_aa[geneID] = current_seq_aa

                        elif strand_by_geneID_dict.get(geneID, "NA") == "-":
                            seq_codons[geneID] = current_seq_codons[::-1]
                            seq_aa[geneID] = current_seq_aa[::-1]
                        
        
                    if strand_by_geneID_dict.get(geneID, "NA") == "+":
                        current_matrix[ix,:] = current_ribo

                    elif strand_by_geneID_dict.get(geneID, "NA") == "-":
                        current_matrix[ix,:] = current_ribo[::-1]
            
            if np.sum(current_matrix) > 0:    
                scikit_mat[geneID] = current_matrix

#        scikit_df = pd.DataFrame(values_list, columns=columns_list)

        return scikit_mat, seq_codons, seq_aa


    def mean_norm(row):
        codon_dens_prof = row.codon_density_profile
        profile_average = np.average(codon_dens_prof)

        return [x/profile_average for x in codon_dens_prof]
    
    #scikit_data_df["mean_norm_codon_density_profile"] = scikit_data_df.apply(mean_norm, axis=1)
    #scikit_data_df["mean_norm_codon_density_profile"] = scikit_data_df['mean_norm_codon_density_profile'].apply(lambda x: x[8:-8])

    strand_by_geneID_dict = gather_strand_by_geneID_dict(genome_gtf)
    scikit_data_dict = import_scikit_data(sralist)
    scikit_data_mat, seq_codons_dict, seq_aa_dict = build_mat_scikit_strandOriented(sralist, scikit_data_dict)

    with open('../data/processed/scikit_mat_pre.pkl', 'wb') as f:
    	pickle.dump(scikit_data_mat, f)

    with open('../data/processed/scikit_codonseq.pkl', 'wb') as f_seq:
        pickle.dump(seq_codons_dict, f_seq)
    

    return scikit_data_mat



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def build_mm_df(sralist):
    """
    dictionary of boolean multi-mapping matrices
    """

    def convert_to_codon(nts_array):
        """
        pysam output is in nucleotides resolution, but scikit_curated_df uses codon resolution.
        This function converts nucleotide arrays to codon length (nts to codon resolution):
        """
        
        nts_array = np.array(nts_array)
        codon_array = np.sum( np.reshape(nts_array, (int(np.floor(len(nts_array)/3)),3) ), 1)/3.

        return codon_array


    def compute_mm(mmdata):
        """
        get per gene average multi-mapping score
        """

        mm_df = pd.DataFrame(columns=['ORF', 'MM'])
        counter = 0

        for gene in mmdata.keys():
            current_matrix = mmdata[gene]
            current_avrg = np.mean( np.sum(current_matrix, 1) / current_matrix.shape[1] )
            mm_df.loc[counter] = [gene, current_avrg]
            counter += 1

        return mm_df


    mm_mat = {}
    mm_pct = {}

    N = len(sralist)

    for ix, dataset in enumerate(sralist):
        samfile = pysam.AlignmentFile(TMP_DIR+'/ambiguous_reads/'+dataset+'_STAR_transcriptome_multi_mapped_sorted.bam', 'rb')
        genes_list = list(samfile.references)
        print(ix, dataset)

        for geneID in genes_list:
            # count the coverage of genomic positions by reads in region.
            # Returns: four array.arrays of the same length in order A C G T
            # The coverage is computed per-base [ACGT]
            cov = samfile.count_coverage(geneID, read_callback='nofilter')
            # Summ all 4 arrays
            cov_sum = np.sum(cov, axis=0)
            #print(geneID, cov_sum)
            codon_cov = convert_to_codon(cov_sum)
            codon_bool = np.asarray([1 if i > 0 else 0 for i in codon_cov])
        
            M = len(codon_bool)

            if ix == 0:
            	mm_mat[geneID] = np.zeros((N,M)) * np.nan
                        
            current_matrix = mm_mat[geneID]
            current_matrix[ix,:] = np.copy(codon_bool)
            mm_mat[geneID] = current_matrix


    mm_avrg = compute_mm(mm_mat)
    #mm_avrg.to_json('yeast_mm.json')
    #mm_avrg.to_csv('yeast_mm.txt', header=True, index=False, sep='\t')

 
    mm_profile = {}
    theta_mm = 5
    for orf in mm_mat.keys():
        current_mat = mm_mat[orf]
        current_bool = np.sum(current_mat, 0) <= theta_mm
        mm_profile[orf] = current_bool

    with open('../data/processed/mm_consensus.pkl', 'wb') as f_mm:
        pickle.dump(mm_profile, f_mm)


    return mm_mat, mm_avrg, mm_profile



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
def filter_normalize_data(ribo_dict, mm_dict, threshold_ribo, threshold_mm):
    """
    final data QC
    """

    final_ribo = {}
    final_mm = {}

    output = pd.DataFrame(columns=['ORF', 'EL', 'MM'])
    counter = 0

    for i in scikit_data.keys():
        if i in mm_data.keys():
            current_ribo = ribo_dict[i]
            current_multi  = mm_dict[i]
            current_el = np.around( np.mean( np.sum(current_ribo, 1) / current_ribo.shape[1]), 3)
            current_mm = np.around( np.mean( np.sum(current_multi, 1) / current_multi.shape[1]), 3)
            output.loc[counter] = [i, current_el, current_mm]
            counter += 1

            if current_el >= threshold_ribo and current_mm <= threshold_mm:
                for j in range(current_ribo.shape[0]):
                    if np.mean(current_ribo[j,:]) > 0:
                        current_ribo[j,:] = current_ribo[j,:] / np.mean(current_ribo[j,:])
                final_ribo[i] = current_ribo
                final_mm[i] = current_multi

    # save output
    output.to_csv('../data/processed/yeast_filter.txt', header=True, index=False, sep='\t')
    pickle.dump(final_mm, open('../data/processed/mm_mat.pkl', 'wb'))
    pickle.dump(final_ribo, open('../data/processed/scikit_mat.pkl', 'wb'))
    







if __name__ == '__main__':


    ## - global variables
    genome 	= '../data/reference/genome-R64-1-1.fa'
    gtf 	= '../data/reference/genes-R64-1-1.gtf' 
    transcriptome = '../data/reference/orf_coding.fasta'		# https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
    rRNA 	= '../data/reference/data/yeast_rRNA.fa'			# non-coding stuff from the yeast gff

    ## - set paths to various directories
    NCBI_DIR = '~/ncbi/public/SRA/'           # set to local NCBI directory
    TMP_DIR = 'TMP/'                        # temporary directory for all intermediate processing files
    BOWTIE_DIR = '../data/indices/bowtie'      # directory for bowtie indices
    STAR_DIR = '../data/indices/STAR'
    SCIKIT_DIR = '../data/indices/scikit'
    STAR_TRANSCRIPTOME_DIR = '../data/indices/transcriptome'


    ## - SET PATH TO LOCAL INSTALL OF MODIFIED SCIKIT-RIBO PACKAGE (JUST EXTRACT DATA, NOT GENERATE PLOTS OR DOWNSTREAM ANALYSES)
    SCIKIT_PATH = 'scikit-ribo-mod/scikit_ribo/'


    ## thresholds for filtering
    theta_ribo = 1    # minimum of 1 read per position on average
    theta_mm   = 0.3  # maximum of 30% multi-mapping
  

    ## - stuff that runs once, like index building
    build_indices(genome, gtf, rRNA, transcriptome)

    
    ## - stuff that runs for each SRA
    # replace this list by glob.glob of the downloaded SRA files in case filenames don't match! 
    list_SRA = ['SRR1002819','SRR1042855','SRR1520311','SRR1688545','SRR1944981','SRR1944982','SRR1944983',
				'SRR2046309','SRR2046310','SRR3493886','SRR3493887','SRR3623557','SRR3623558','SRR3945926',
				'SRR3945928','SRR4000288','SRR4000289','SRR5008134','SRR5008135','SRR5090936']


    for i in range( len(list_SRA) ):
        
        current_sra = list_SRA[i]

        download_SRA(current_sra)				# 
        trim_adapters(current_sra)				# this runs! 
        map_reads(current_sra)					# this runs!
        run_scikit_ribo(current_sra, genome, gtf)			
        run_multimapping(current_sra)


    ## - stuff that runs once, AFTER the SRA files have been processed
    ## - generate scikit df(s)
    scikit_data = parse_dataframes(gtf, list_SRA)
 
    ## - generate multimapping df(s)
    mm_data, mm_avrg, mm_consensus = build_mm_df(list_SRA)

    ## - filter/curate final list of genes
    filter_normalize_data(scikit_data, mm_data, theta_ribo, theta_mm)

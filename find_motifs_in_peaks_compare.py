#!/usr/bin/env python
# coding: utf-8

# In[1]:
import pyBigWig
import numpy as np
import pandas as pd

# In[2]:
def shap_threshold(statfile, quantile='99%'):
    score_stats = pd.read_csv(statfile, sep='\t', index_col=0, header=None)
    threshold=score_stats.loc[quantile].values
    return threshold

def peak(bed, n):
    pid = bed.index[n]
    chrom = bed.chr.iloc[n]
    start = int(bed.start.iloc[n])
    end = int(bed.end.iloc[n])
    return pid, chrom, start, end

def overlap(start1, end1, start2, end2):
    return min(end1,end2) >= max(start1,start2) and max(end1,end2) >= min(start1,start2)

def overlap_mat(motifs):
    n = len(motifs)
    ranges = [(m[1],m[2]) for m in motifs]
    o_mat = np.zeros((n,n))
    
    from itertools import product
    for i,j in product(range(n),range(n)):
        o_mat[i,j] = overlap(ranges[i][0], ranges[i][1], ranges[j][0],ranges[j][1])
    return o_mat

def relm_score(motif):
    return motif[4]*max([motif[5],motif[6]])

def filter_motifs(motifs):
    n = len(motifs)
    relm_scores = [relm_score(m) for m in motifs]
    o_mat = overlap_mat(motifs)
    filter_list=list(range(len(motifs)))
    sorted_idx=np.argsort(relm_scores)[::-1]
    for m in range(n):
        if filter_list[sorted_idx[m]]==sorted_idx[m]:
            overlaps = [i for i, x in enumerate(o_mat[sorted_idx[m]]) if x]
            for o in overlaps:
                filter_list[o]=sorted_idx[m]
    return list(np.unique(filter_list))

def process_peak(args):
    p, sig_peaks, motif_bb, shap1, shap2, threshold1, threshold2 = args
    
    pid, chrom, start, end = peak(sig_peaks,p)
    print('Peak {}: {}:{}-{}'.format(pid,chrom,start,end))

    peak_motifs = motif_bb.entries(chrom, start, end)

    threshold_motifs = []
    if peak_motifs:
        for motif in peak_motifs:
            m_start, m_end, m_type, match_score = motif[0], motif[1], motif[2].split('\t')[0], float(motif[2].split('\t')[7])
            m_score_1 = float(sum(shap>threshold1 for shap in shap1.values(chrom, m_start, m_end))/(m_end-m_start))
            m_score_2 = float(sum(shap>threshold2 for shap in shap2.values(chrom, m_start, m_end))/(m_end-m_start))
            if any([m_score_1>0.25, m_score_2>0.25]):
                threshold_motifs += [[chrom, m_start, m_end, m_type, match_score, m_score_1, m_score_2]]

    filtered_motifs = [threshold_motifs[x]+[pid] for x in filter_motifs(threshold_motifs)]
    if filtered_motifs:
        print('Found motifs: {}'.format([m[3] for m in filtered_motifs]))
        return filtered_motifs
    
def main(sig_peaks, bw_1, bw_2, stats1, stats2, genome):
    if genome == 'mm10':
        motif_bb = pyBigWig.open("https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/mm10.archetype_motifs.v1.0.bb")
    elif species == 'hg38':
        motif_bb = pyBigWig.open("https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bb")
    else:
        raise NameError('Only mm10 and hg38 are supported at this time')
        
    shap1 = pyBigWig.open(bw_1)
    shap2 = pyBigWig.open(bw_2)

    threshold1 = shap_threshold(stats1, quantile='99%')
    threshold2 = shap_threshold(stats2, quantile='99%')
    #threshold1 = 0.01
    #threshold2 = 0.01

    # In[5]:
    #sig_peaks.columns = ['peak_id', 'chr', 'start', 'end', 'Gene 1', 'Gene 1 dist', 'Gene 2', 'Gene 2 dist',
    #       'Gene 1 FC', 'Gene 2 FC', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat',
    #       'pvalue', 'padj', 'maxDeltaFC', 'maxDeltaGene', 'maxDeltaDist',
    #       'Direction', 'CRE half sites "TGACG and CGTCA"', 'Full CRE "TGACGTCA"',
    #       'AP1 "TGAGTCA and TGACTCA"', 'Ebox "CAGCTG"']

    sig_peaks = sig_peaks.set_index('peak_id')
    #sig_peaks.index=[int(x[4:]) for x in sig_peaks.index]
    sig_peaks=sig_peaks.sort_index()


    # In[ ]:
    data = [(i, sig_peaks, motif_bb, shap1, shap2, threshold1, threshold2) for i in range(len(sig_peaks))]
    filtered_motifs=[process_peak(d) for d in data]


    # In[101]:
    motif_output = [m for peak in filtered_motifs if peak for m in peak ]
    motif_df = pd.DataFrame(motif_output)
    motif_df.columns=['chr','start','end','motif_cluster','match_score','shap1','shap2','peak_id']

    # In[163]:
    motif_df.to_csv('filtered_BPNet_motifs.bed', sep='\t',index=False)

    # In[164]:
    peak_dummies = pd.get_dummies(motif_df,columns=['motif_cluster'], prefix='', prefix_sep='').groupby('peak_id').sum().iloc[:,3:]
    peak_dummies['n_motifs']=peak_dummies.iloc[:,2:].sum(axis=1)

    peaks_output = sig_peaks.join(peak_dummies)
    peaks_output.iloc[:,-peak_dummies.shape[1]+2:]=peaks_output.iloc[:,-peak_dummies.shape[1]+2:].fillna(0)
    peaks_output['peak_id']=sig_peaks.index

    # In[165]:
    peaks_output.to_csv('20210622_youngE2_sigpeaks_BPNet_motifs_99.csv', index=False)


# In[3]:
if __name__ == '__main__':

    #aging
    #sig_peaks=pd.read_csv('/Users/will/Documents/basepairmodels/scg/sigpeaks/20210615_aging_sigpeaks_gene_associations.csv')
    #bw1 = '/Users/will/Documents/basepairmodels/scg/young_2h/interpretations/output/veh_v2a2_2_score.bw'
    #bw2 = '/Users/will/Documents/basepairmodels/scg/aged_2h_v2/interpretations/output/veh_v2a4_score.bw'

    #youngE2
    sig_peaks=pd.read_csv('/Users/will/Documents/basepairmodels/scg/sigpeaks/20210615_youngE2_sigpeaks_gene_associations.csv')
    bw1 = '/Users/will/Documents/basepairmodels/scg/young_2h/interpretations/output/veh_v2a2_2_score.bw'
    bw2 = '/Users/will/Documents/basepairmodels/scg/young_2h/interpretations/output/E2_v2a2_2_score.bw'

    #agedE2
    #sig_peaks=pd.read_csv('/Users/will/Documents/basepairmodels/scg/sigpeaks/20210622_agedE2_sigpeaks_gene_associations.csv')
    #bw1 = '/Users/will/Documents/basepairmodels/scg/aged_2h_v2/interpretations/output/veh_v2a4_score.bw'
    #bw2 = '/Users/will/Documents/basepairmodels/scg/aged_2h_v2/interpretations/output/E2_v2a4_score.bw'

    stats1 = bw1[:-3]+'.stats.txt'
    stats2 = bw2[:-3]+'.stats.txt'

    genome = 'mm10'

    main(sig_peaks, bw1, bw2, stats1, stats2, genome)
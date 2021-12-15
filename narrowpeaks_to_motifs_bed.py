import pyBigWig
import numpy as np
import pandas as pd
from tqdm import tqdm

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
    return motif[4]*motif[5]

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
    p, peaks, motif_bb, shap, threshold = args
    
    pid, chrom, start, end = peak(peaks,p)
    #print('Peak {}: {}:{}-{}'.format(pid,chrom,start,end))

    peak_motifs = motif_bb.entries(chrom, start, end)

    threshold_motifs = []
    if peak_motifs:
        for motif in peak_motifs:
            m_start, m_end, m_type, match_score = motif[0], motif[1], motif[2].split('\t')[0], float(motif[2].split('\t')[7])
            m_score = float(sum(shap>threshold for shap in shap.values(chrom, m_start, m_end))/(m_end-m_start))
            if m_score>0.25:
                threshold_motifs += [[chrom, m_start, m_end, m_type, match_score, m_score]]

    filtered_motifs = [threshold_motifs[x]+[pid] for x in filter_motifs(threshold_motifs)]
    if filtered_motifs:
        #print('Found motifs: {}'.format([m[3] for m in filtered_motifs]))
        return filtered_motifs

def number(str):
    r=''
    for i in str:
        if i.isdigit():
            r+=i
    return int(r)

def main(peaksfile, shap_bw, threshold, genome, motifs_outfile, peaks_outfile):
    if genome == 'mm10':
        motif_bb = pyBigWig.open("https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/mm10.archetype_motifs.v1.0.bb")
    elif species == 'hg38':
        motif_bb = pyBigWig.open("https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bb")
    else:
        raise NameError('Only mm10 and hg38 are supported at this time')

    narrowpeaks = pd.read_csv(peaksfile, sep='\t', header=None)
    narrowpeaks.columns = ['chr', 'start', 'end', 'summit_id', 'score', 'strand', 'signal', 'p', 'q', 'peak']
    narrowpeaks['peak_id'] = [number(id.split('_')[2]) for id in narrowpeaks['summit_id']]

    peak_coords = narrowpeaks.groupby('peak_id').first()[['chr','start','end']]

    shap = pyBigWig.open(shap_bw)

    data = [(i, peak_coords, motif_bb, shap, threshold) for i in range(len(peak_coords))]
    filtered_motifs=[process_peak(d) for d in tqdm(data)]

    motif_output = [m for peak in filtered_motifs if peak for m in peak ]
    motif_df = pd.DataFrame(motif_output)
    motif_df.columns=['chr','start','end','motif_cluster','match_score','shap','peak_id']
    motif_df.to_csv(motifs_outfile, sep='\t',index=False)

    peak_dummies = pd.get_dummies(motif_df,columns=['motif_cluster'], prefix='', prefix_sep='').groupby('peak_id').sum().iloc[:,3:]
    peak_dummies['n_motifs']=peak_dummies.iloc[:,1:].sum(axis=1)

    peaks_output = peak_coords.join(peak_dummies)
    peaks_output.iloc[:,-peak_dummies.shape[1]+1:]=peaks_output.iloc[:,-peak_dummies.shape[1]+1:].fillna(0)
    peaks_output['peak_id']=peak_coords.index
    
    peaks_output.to_csv(peaks_outfile, index=False)

if __name__ == '__main__':

    peaksfile = '/Users/will/Documents/basepairmodels/scg/aged_2h_v2/data/aged2h.narrowpeak'
    shap_bw = '/Users/will/Documents/basepairmodels/scg/aged_2h_v2/interpretations/output/E2_v2a4_score.bw'
    motifs_outfile = '/Users/will/Documents/basepairmodels/scg/aged_2h_v2/interpretations/output/E2_v2a4_99_motifs.bed'
    peaks_outfile = '/Users/will/Documents/basepairmodels/scg/aged_2h_v2/interpretations/output/E2_v2a4_99_peak_motifs.csv'

    genome = 'mm10'
    
    if 1:
        shap_stats = shap_bw[:-3]+'.stats.txt'
        quantile = '99%'
        threshold = shap_threshold(shap_stats, quantile)
    else:
        threshold = 0.01#shap_threshold(shap_stats, quantile)

    main(peaksfile, shap_bw, threshold, genome, motifs_outfile, peaks_outfile)
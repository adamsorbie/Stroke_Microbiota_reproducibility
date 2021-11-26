#!/usr/bin/env python3
import pandas as pd
import numpy as np 

def normalize(asvtab):
    sum_row = asvtab.sum()
    rel_abund = asvtab * 100 / sum_row
    return rel_abund

def filter_abund(asvtab, abun_thresh=0.25, prev_thresh=1):
    asvtab_filt = asvtab[(asvtab > abun_thresh).sum(axis=1) >= prev_thresh]
    asv_idx = asvtab_filt.index.tolist()
    return(asv_idx)

def threshold_asvtab(asvtab):
    asvtab_rel = normalize(asvtab)
    filt_asvs = filter_abund(asvtab_rel)
    asvtab_out = asvtab[asvtab.index.isin(filt_asvs)]
    return(asvtab_out, filt_asvs)

asvtab = pd.read_csv("feature-table.tsv")


asvtab_filt, asv_filt_list = threshold_asvtab(asvtab)


asvtab_filt.to_csv("feature-table-filt.tsv", sep="\t")

asv_filt_list_out = pd.DataFrame(asv_filt_list)
asv_filt_list_out.to_csv('keep_asvs.tsv', index=False, sep="\t")





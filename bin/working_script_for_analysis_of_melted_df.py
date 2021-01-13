#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 01:21:59 2020

@author: ptruong
"""

os.chdir("/home/ptruong/git/bayesMS/bin")

import os 

import pandas as pd
import numpy as np 

from read_triqler_output import read_triqler_protein_output_to_df
from triqler_output_df_melter import melt_spectronaut_triqler_formatted, melt_triqler_output


spec = pd.read_pickle(r'spectronaut.pkl')
spec = spec.rename(columns={'S03:S04_R05': 'S03:S03_R04'})
triq = pd.read_pickle(r'triqler.pkl')

triq = read_triqler_protein_output_to_df('proteins.3vs8.tsv')


triq = melt_triqler_output(triq)
spec = melt_spectronaut_triqler_formatted(spec)


samples = ["S0"+str(i) for i in range(1,10)] + ["S10"]


triq.columns

triq[triq["protein_id_posterior_error_prob"] > 0.01] #Both triqler and spec are already 1% FDR tresholded.


# Identified protein through all samples
print("%s : %i" % ("triqler number of protein ids", int(sum(triq.protein_id_posterior_error_prob < 0.01))))
print("%s : %i" % ("spectronaut number of protein ids", int(len(spec.dropna()))))

# Identified protein for species across all samples
triq_count = len(triq[triq.specie == "ARATH"])
spec_count = len(spec[spec.specie == "ARATH"].dropna())

print("%s : %i" % ("triqler ARATH number of protein ids", triq_count))
print("%s : %i" % ("spectronaut ARATH number of protein ids", spec_count))

triq_count = len(triq[triq.specie == "CAEEL"])
spec_count = len(spec[spec.specie == "CAEEL"].dropna())

print("%s : %i" % ("triqler CAEEL number of protein ids", triq_count))
print("%s : %i" % ("spectronaut CAEEL number of protein ids", spec_count))

triq_count = len(triq[triq.specie == "HUMAN"])
spec_count = len(spec[spec.specie == "HUMAN"].dropna())

print("%s : %i" % ("triqler HUMAN number of protein ids", triq_count))
print("%s : %i" % ("spectronaut HUMAN number of protein ids", spec_count))

# Create an identification matrix (df) for each all samples, that can have species as input

samples

specie = None
def count_melted_for_all_samples(df, specie = None):
    """
    df = melted triq or spec.
    """
    counts = []
    for i in samples:
        count_df = df[df["sample"] == i]
        if specie == "ARATH":
            count_df = count_df[count_df["specie"] == "ARATH"].dropna()
        elif specie == "HUMAN":
            count_df = count_df[count_df["specie"] == "HUMAN"].dropna()
        elif specie == "CAEEL":
            count_df = count_df[count_df["specie"] == "CAEEL"].dropna()
        count = len(count_df.dropna())
        counts.append(count)
    return counts

count_melted_for_all_samples(triq, specie = None)
count_melted_for_all_samples(triq, specie = "HUMAN")
count_melted_for_all_samples(triq, specie = "CAEEL")
count_melted_for_all_samples(triq, specie = "ARATH")

count_melted_for_all_samples(spec, specie = None)
count_melted_for_all_samples(spec, specie = "HUMAN")
count_melted_for_all_samples(spec, specie = "CAEEL")
count_melted_for_all_samples(spec, specie = "ARATH")



# Normalization

from normalize_melted import normalize_within_sample, get_ratios_from_normalized_melted_df

df_spec = normalize_within_sample(spec)
df_triq  = normalize_within_sample(triq)

# Explore normalization


get_ratios_from_normalized_melted_df(df_triq, "HUMAN")
get_ratios_from_normalized_melted_df(df_triq, "ARATH")
get_ratios_from_normalized_melted_df(df_triq, "CAEEL")



# Compute True FC on samples


def get_log2FC_ratio_matrix(specie_array):
    """
    Inpur specie array with mixture ratios, and generate log2FC matrix.
    """
    samples = ["S0"+str(i) for i in range(1,10)] + ["S10"] 
    FC_ratios = []
    for i in range(10):
        sample_1 = specie_array[i]
        FC_ratio = []
        for j in range(10):
            sample_2 = specie_array[j]
            FC_ratio.append(np.log2(sample_1) - np.log2(sample_2))
        FC_ratios.append(FC_ratio)
    df_FC_ratios = pd.DataFrame(FC_ratios, index = samples, columns = samples)
    return df_FC_ratios

def get_log2FC_ratio_matrices():
    ARATH = np.array([0.5,  0.5,  0.5,  0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    CAEEL = np.array([ 0.5,  0.25,  0.125,  0.0625,  0.031,  0.0155, 0.008,  0.004,  0.002,  0.000001])
    HUMAN = np.array([ 00.000001,  0.25, 0.375,  0.4375,  0.469,  0.4845,  0.492, 0.496, 0.498,  0.5])
    
    ARATH_FC_matrix = get_log2FC_ratio_matrix(ARATH)
    CAEEL_FC_matrix = get_log2FC_ratio_matrix(CAEEL)
    HUMAN_FC_matrix = get_log2FC_ratio_matrix(HUMAN)
    return ARATH_FC_matrix, CAEEL_FC_matrix, HUMAN_FC_matrix


ARATH_FC_matrix, CAEEL_FC_matrix, HUMAN_FC_matrix = get_log2FC_ratio_matrices()
    
# Computing FC diff-exp for triq and spec
samples = ["S0"+str(i) for i in range(1,10)] + ["S10"] 

df = triq #variable
df = spec #variable

sample = samples[0] #Variable
specie = "HUMAN" #variable

df_sample = df[df["sample"] == sample]
df_sample_specie = df_sample[df_sample["specie"] == specie]

proteins = df_sample_specie.protein.unique()

protein = proteins[0] #iteration variable (if the protein is diff exp)

df_sample_specie_protein = df_sample_specie[df_sample_specie["protein"] == protein].value






# Compute similar thing for proteinXvsY.csv - the posterior diff exp.








































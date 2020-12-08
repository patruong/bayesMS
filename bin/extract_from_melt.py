#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 00:55:45 2020

@author: ptruong
"""




len(spec.protein.unique())
protein_list = spec.protein.unique()

protein = protein_list[0]

protein_df = spec[spec.protein == protein]
protein_df[protein_df["sample"] == "S01"]



#Check normalization of melted

len(triq.protein.unique())
protein_list = triq.protein.unique()
protein = protein_list[0]

protein_df = triq[triq.protein == protein].drop("peptide", axis = 1)
protein_run_df = protein_df[protein_df.run == "R01"]
protein_run_df.value


## first lets try to get matthews code working for the plotting.

















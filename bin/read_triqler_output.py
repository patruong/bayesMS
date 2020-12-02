#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 14:43:35 2020

@author: ptruong
"""

import pandas as pd
import numpy as np

filename = "proteins.1vs2.tsv"
f = open(filename, "r")

header = f.readline().split("\n")[0].split("\t")

len_header = len(header)

lines = []
for line in f:
    line = line.split("\n")[0].split("\t")
    vals = line[:len_header - 1]
    peptides = line[(len_header - 1 ):]
    peptides = ";".join(peptides)
    vals.append(peptides)
    lines.append(vals)

df = pd.DataFrame(np.array(lines).T)




#!/usr/bin/env python3

import os
import pandas as pd

summaries = []

for root, _, files in os.walk("."):
    for file in files:
        if file.endswith('.tsv'):
            summ_f = os.path.join(root,file)
            summ = pd.read_table(summ_f)
            summaries.append(summ)

summ_m = pd.concat(summaries)
summ_m.to_csv("classification_summary.tsv", sep="\t",index=False)


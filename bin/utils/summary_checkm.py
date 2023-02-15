#!/usr/bin/env python3

import os
import pandas as pd

contents = []

for root, _, files in os.walk("."):
    for f in files:
        if f == "checkm_results.tsv":
            path = os.path.join(root, f)
            df = pd.read_table(path)

            contents.append(df)

summary = pd.concat(contents)
print(summary)
summary.to_csv("summary_checkm.tsv", sep="\t", index=False)

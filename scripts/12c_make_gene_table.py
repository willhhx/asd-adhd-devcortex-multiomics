#!/usr/bin/env python3
"""Gene TSS table (gene, chrom, tss; hg19) from UCSC refGene.txt.gz for the nearest-TSS fallback of step 12.
Download: curl -L -o data/raw/refGene.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
"""
import argparse

import pandas as pd

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--refgene", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
d = pd.read_csv(a.refgene, sep="\t", header=None, usecols=[2, 3, 4, 5, 12], names=["chrom", "strand", "s", "e", "gene"])
d = d[d["chrom"].str.match(r"^chr([0-9]+|X)$")]
d["chrom"] = d["chrom"].str.replace("chr", "", regex=False)
d["tss"] = d["s"].where(d["strand"] == "+", d["e"])
d[["gene", "chrom", "tss"]].drop_duplicates().to_csv(a.out, sep="\t", index=False)
print(f"{d['gene'].nunique()} genes, {len(d)} TSS rows -> {a.out}")

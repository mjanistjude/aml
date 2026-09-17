import pandas as pd
import sys
import os

df=pd.read_csv(sys.argv[1],sep="\t")
df["chr2"] = df["ALT"].str.extract(r"chr(.*?):")
df["pos2"] = df["ALT"].str.extract(r":(\d+)")
df["str1"] = df["INFO"].str.extract(r"STRANDS=(.)", expand=False).map({"+":"0","-":"1"})
df["str2"] = df["INFO"].str.extract(r"STRANDS=.(.)", expand=False).map({"+":"0","-":"1"})
df["chr1"] = df["#CHROM"].str.replace("chr","")
df=df[["chr1","POS","str1","chr2","pos2","str2"]]
df.columns=["chr1","pos1","str1","chr2","pos2","str2"]
outfile=os.path.basename(sys.argv[1]).split("_extracted.tsv")[0]
df.to_csv(sys.argv[2]+outfile+"_input.tsv",sep="\t",index=False)

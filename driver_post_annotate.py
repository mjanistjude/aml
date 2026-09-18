#source activate agfusion
import orientation_check as oc
import prep_agfusion as agf
import pandas as pd
import frame
import sys
import os

infile=sys.argv[1]
outdir=sys.argv[2]
name=sys.argv[3]
hgref=sys.argv[4]
mittledb=sys.argv[5]

outfile=outdir+"/"+name+".agfusinput.tsv"
orifile=outdir+"/"+name+".ag.tsv"
mdf=pd.read_csv(mittledb,sep="\t")["Fusion Gene"].to_list()
df=agf.driver(outfile,infile,outdir,name)
print ("agfusion: done")
df,hgdf=oc.driver(hgref,orifile)
print ("orientation_done")
df["status"]=df.apply(lambda x: frame.driver_frame(hgdf,hgref,x["Gene_St"],x["transcript1"],x["Gene_End"],x["transcript2"],int(x["POS"]),int(x["endpos"])), axis=1)
print ("writing to file")
df["mittleman"] = df.apply(lambda x: "Y" if pd.notna(x["Gene_St"]) and pd.notna(x["Gene_End"]) and (str(x["Gene_St"]) + "::" + str(x["Gene_End"]) in mdf) else "-",axis=1)
df.to_csv(outdir+"/"+name+"_annotated.sag.tsv",sep="\t",index=False)

df=df[["#CHROM","POS","endchr","endpos","SV_length","BAR","splitAsplitB","Gene_St","site1","Gene_End","site2","Fusion_effect","coverage_bp1","coverage_bp2","ratioA","ratioB","mittleman","GroupID"]]
df.columns=["chrA","posA","chrB","posB","size","BAR","sf","GeneA","exon/intronA","GeneB","exon/intronB","Frame","Total Reads A","Total Reads B","RatioA","RatioB","Tags","GroupID"]

df[["sfa", "sfb"]] = df["sf"].str.split(",",expand=True)
df["Tags"] = df["Tags"].replace("Y", "Mitelman")
df.drop(columns=["sf"], inplace=True)

df.to_csv(outdir+"/"+name+"_aperture_review.tsv",sep="\t",index=False)


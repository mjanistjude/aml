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
df.to_csv(outdir+"/"+name+".sag.tsv",sep="\t",index=False)


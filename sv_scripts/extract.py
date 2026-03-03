import pandas as pd
import glob
import sys

#files=sorted(glob.glob("*.vcf"))

def bed(dfref,svbed):
	cols=["gene","chr","start","end"]
	df=pd.read_csv(svbed,sep="\t",names=cols)
	df=df[["chr","start","end"]]
	arr=[]
	for rec in dfref.to_dict("records"):
		df1=df.loc[ (df["chr"]==rec["#CHROM"]) & (df["start"]<=int(rec["POS"])) & (df["end"]>=int(rec["POS"])) ]
		if len(df1)>0:
			rec["bedst"]="Yes"
		else:
			rec["bedst"]="No"
		df1=df.loc[ (df["chr"]==rec["endchr"]) & (df["start"]<=int(rec["end"])) & (df["end"]>=int(rec["end"])) ]
		if len(df1)>0:
			rec["bedend"]="Yes"
		else:
			rec["bedend"]="No"
		arr.append(rec)
	dfref=pd.DataFrame(arr)
	dfref=dfref.loc[ (dfref["bedst"]=="Yes") | (dfref["bedend"]=="Yes") ]
	return dfref

def extract_info_field(info, field):
	for item in info.split(';'):
		if item.startswith(field + '='):
			return item.split('=')[1]
	return -1

def exp_vcf(df):
	for c in ["REFQUA","VARQUA","REFKMER","VARKMER","BPSEQQUA"]:
		df[c] = df['INFO'].apply(lambda x: extract_info_field(x, c))
		df[c] = pd.to_numeric(df[c], errors='coerce')
	sample_columns = ['GT', 'SR', 'PE', 'REFSR', 'VARSR', 'BAR', 'UBAR']
	df[sample_columns]=df['SAMPLE'].str.split(':', expand=True)
	for c in sample_columns[1:]:
		df[c] = pd.to_numeric(df[c])
	return df


def endchr(x):
	return "chr"+x["ALT"].split(":")[0].split("chr")[1]

def end(x):
	return int(x["ALT"].split(":")[1].split("[")[0].split("]")[0])

def getval(x,col):
	val=x.split(";")
	val=[v for v in val if col in v]
	if len(val)==0:
		if col=="HOMLEN":
			return 0
		return ""
	val=val[0].split("=")[1]
	if col in ["STRANDS","HOMSEQ"]:
		return val
	return float(val)

def filter_vcf(df):
	def revlab(x):
		if x["BAR"]>2 and x["SR"]>10 and x["VARSR"]>10 and x["VARQUA"]>3 and x["REFQUA"]>3 and x["VARKMER"]>100 and x["BAR"]>x["UBAR"]:
			return ""
		else:
			return "Need Review"
	cols=list(df)[:-1]+["FNUMS"]
	df.columns=cols
	for col,num in zip(["SR","PE","REFSR","VARSR","BAR","UBAR"],range(1,7)):
		df[col]=df.apply(lambda x: float(x["FNUMS"].split(":")[num]), axis=1)
	for col in ["REFQUA","VARQUA","REFKMER","VARKMER","BPSEQQUA","STRANDS","HOMLEN","HOMSEQ"]:
		df[col]=df.apply(lambda x: getval(x["INFO"],col), axis=1)
	df["Review"]=df.apply(lambda x: revlab(x),axis=1)
	return df


files=[sys.argv[1]]
svbed=sys.argv[2]
for ff in files:
	sample=ff.split(".")[0]
	header=""
	f=open(files[0],"r")
	for line in f:
		if line.startswith("#"):
                	header+=line
	                last=line.strip().split("\t")

	df=pd.read_csv(ff,comment="#",sep="\t",names=last)
	df=df.loc[df["FILTER"]=="PASS"]
	df.to_csv(sample+".pass.tsv",sep="\t",index=False)
	print (df)
	df["endchr"]=df.apply(lambda x: endchr(x),axis=1)
	df["end"]=df.apply(lambda x: end(x),axis=1)
	df=bed(df,svbed)
	df = df.drop(columns=['endchr','end','bedst','bedend'])
	out_f=open(ff.split(".")[0]+"_extracted.vcf","w")
	out_f.write(header)
	df.to_csv(out_f, sep="\t", mode='a', index=False, header=False)
	f.close()
	out_f.close()
	df=filter_vcf(df)
	print (list(df))
	df=df[["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","FNUMS","Review"]]
	print (df)
	print (sys.argv[3])
	df.to_csv(sys.argv[3], sep="\t", index=False)


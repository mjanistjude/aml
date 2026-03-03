import pandas as pd
import glob
import os

def parse(df,outname):
	df=df[["Gene_St","ID","GroupID","SV_type","Gene_End","site1","site2","GA_GB","REF","POS","ALT","endpos","INFO"]]
	df=df.loc[~((df["site1"].str.contains("IGR:"))|(df["site2"].str.contains("IGR:")))]
	df=df.loc[~(df["Gene_St"].isna() | df["Gene_End"].isna())]
	#df=df.loc[df["GA_GB"]=="Y"]
	df.to_csv(outname,sep="\t",index=False,header=None)
	return

def run(infile,outdir):
	cmd="agfusion batch -f %s -a bellerophontes -db /clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/SV_output/agfusion/agfusion.homo_sapiens.75.db -o %s"%(infile,outdir)
	print (cmd)
	os.system(cmd)
	return

def create_outfile(outdir):
	dirs=[entry.name for entry in os.scandir(outdir) if entry.is_dir()]
	alldfs=[]
	if not dirs:
		print(f"Warning: No subdirectories found in {outdir}.")
		return pd.DataFrame()
	#print(f"--- Found {len(dirs)} sample directories to process. ---")
	for d in dirs:
		file_paths = glob.glob(os.path.join(outdir, d, "*.fusion_transcripts.csv"))
		if file_paths:
			df = pd.read_csv(file_paths[0])
			df['fus_ID'] = d
			df["Gene_St"]=d.split("_")[0].split("-")[0]
			df["POS"]=int(d.split("_")[0].split("-")[-1])
			df["Gene_End"]=d.split("_")[1].split("-")[0]
			df["endpos"]=int(d.split("_")[1].split("-")[-1])
			alldfs.append(df)
		else:
			print(f"Warning: No matching file found for sample {d}. Skipping.")
	if alldfs:
		combined_df = pd.concat(alldfs, ignore_index=True)
		return combined_df
	else:
		print("\nError: No dataframes were successfully loaded.")
		return pd.DataFrame()

def driver(agfile,annfile,outdir,name):
	df=pd.read_csv(annfile,sep="\t")
	parse(df,agfile)
	run(agfile,outdir)
	combdf=create_outfile(outdir)
	#print (df[["Gene_St","POS","Gene_End","endpos"]])
	df=df.merge(combdf,on=["Gene_St","POS","Gene_End","endpos"],how="left")
	df.to_csv(outdir+"/"+name+".ag.tsv",sep="\t",index=False)
	print (df)
	#print (df[["Gene_St","POS","Gene_End","endpos","Fusion_effect"]])
	return df

if __name__ == "__main__":
	outdir="TB-21-45626_S1"
	annfile=outdir+"_annotated.group.tsv"
	agfile=outdir+".agfusinput.tsv"
	driver(agfile,annfile,outdir)
	files=glob.glob("*.group.tsv")
	for ff in files:
		annfile=ff
		outdir=ff.split("_annotated.group.tsv")[0]
		agfile=outdir+".agfusinput.tsv"
		driver(agfile,annfile,outdir)
		print (ff)

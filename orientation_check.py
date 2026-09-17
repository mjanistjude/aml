import pandas as pd

def read_gtf(path):
    df=pd.read_csv(path,sep="\t",names=["chr","source","feature","start","end","score","strand","frame","attribute"])
    df["gene_id"] = df["attribute"].str.extract(r'gene_id\s+"([^"]+)"', expand=False)
    df["transcript_id"] = df["attribute"].str.extract(r'transcript_id\s+"([^"]+)"', expand=False)
    df["gene_name"] = df["attribute"].str.extract(r'gene_name\s+"([^"]+)"', expand=False)
    df["exon_num"] = df["attribute"].str.extract(r'exon_number\s+"([^"]+)"', expand=False)
    #print (df)
    return df

def check_fusion_sense(gene_a_strand, gene_b_strand, vcf_alt_string):
    """
    Returns True if fusion is SENSE (canonical), False if ANTISENSE.
    gene_a_strand: '+' or '-'
    gene_b_strand: '+' or '-'
    vcf_alt_string: The string from the ALT column (e.g. "T[chr2:1000[")
    """
    if gene_a_strand not in ["-","+"] or gene_b_strand not in ["-","+"]:
        return "IGR/Gene not found"
    # 1. Determine which side of Gene A (REF) is kept
    # If the bracket is AFTER the base (e.g. T[...), we keep the LEFT side of Gene A
    keep_a_left = vcf_alt_string[0].isalpha() or (vcf_alt_string[0] not in ['[', ']'])
    
    # 2. Determine which side of Gene B (ALT) is kept
    # '[' means we keep the RIGHT side of B. ']' means we keep the LEFT side of B.
    # Note: We just look for the bracket character present in the string.
    if '[' in vcf_alt_string:
        keep_b_right = True
    else: # ']' is in string
        keep_b_right = False

    # 3. Check logic for Gene A (Must keep 5' end)
    # If A is (+), 5' is Left. If A is (-), 5' is Right.
    valid_a = (gene_a_strand == '+' and keep_a_left) or \
              (gene_a_strand == '-' and not keep_a_left)

    # 4. Check logic for Gene B (Must keep 3' end)
    # If B is (+), 3' is Right. If B is (-), 3' is Left.
    valid_b = (gene_b_strand == '+' and keep_b_right) or \
              (gene_b_strand == '-' and not keep_b_right)

    if valid_a and valid_b:
        return "SENSE (Canonical)"
    else:
        return "ANTISENSE (Non-canonical)"

def ref_strand(gene,site,trans,pos,df):
	if site is None or (isinstance(site, float)):
		return ""
	if "IGR" in site:
		return "IGR"
	df_filtered=df.loc[df["gene_name"]==gene]
	if df_filtered.empty:
		df_filtered=df.loc[df["transcript_id"]==trans]
		if df_filtered.empty:
			#print (gene)
			df=df.loc[(df["start"]<=pos) & (df["end"]>=pos)]
			if df.empty:
				return "Gene not found"
			strand=list(set(df["strand"].to_list()))
			if len(strand)==1:
				return strand[0]
			else:
				return "Gene not found"
		#return df_filtered["strand"].to_list()[0]
	strand=df_filtered["strand"].to_list()[0]
	return strand

def orient(df,hgdf):
	df["ref_gene_start_strand"]=df.apply(lambda x: ref_strand(x["Gene_St"],x["site1"],x["transcript1"],int(x["POS"]),hgdf), axis=1)
	df["ref_gene_end_strand"]=df.apply(lambda x: ref_strand(x["Gene_End"],x["site2"],x["transcript2"],int(x["endpos"]),hgdf), axis=1)
	df["orientation"]=df.apply(lambda x: check_fusion_sense(x["ref_gene_start_strand"], x["ref_gene_end_strand"], x["ALT"]), axis=1)
	print (df)
	df.to_csv("orient.tsv",sep="\t",index=False)
	return df

def driver(hgref,infile):
	hgdf=read_gtf(hgref)
	df=pd.read_csv(infile,sep="\t")
	print (list(df))
	df=orient(df,hgdf)
	return df,hgdf

if __name__ == "__main__":
	driver("hg19.refGene.gtf","TB-20-06187_S2.ag.tsv")



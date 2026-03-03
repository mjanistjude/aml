import pandas as pd
import splitreads1
import numpy as np
import sys
import re

df=pd.read_csv(sys.argv[1],sep="\t")

def getloc(annotation_string):
    """
    Classifies the annotation string, returning Exon number, IGR with gene name,
    UTR types with gene name, or the calculated Intron number.
    """
    if annotation_string is None:
        return "-"
    if not isinstance(annotation_string, str):
        annotation_string = str(annotation_string)
    # Regex to extract gene name from UTR/Intron strings:
    # Looks for 'of ' followed by the gene name, followed by '(', '+', or '-'
    UTR_INTRON_GENE_PATTERN = r"of\s+([\w\d]+)(?:\(|\+|-)"

    # --- 1. Highest Priority: Exon X ---
    match_exon = re.match(r"Exon (\d+)", annotation_string, re.IGNORECASE)
    if match_exon:
        return f"Exon {match_exon.group(1)}"


    # --- 2. UTR and Intron Gene Name Extraction Helper ---
    # This block handles gene name extraction for 5'-UTR, 3'-UTR, and Intron
    gene_name = "Unknown Gene"
    match_gene = re.search(UTR_INTRON_GENE_PATTERN, annotation_string, re.IGNORECASE)
    if match_gene:
        gene_name = match_gene.group(1)


    # --- 3. Second Priority: UTR Annotation Types (now includes gene name) ---
    if "5'-UTR" in annotation_string:
        return f"5'-UTR ({gene_name})"
    
    if "3'-UTR" in annotation_string:
        return f"3'-UTR ({gene_name})"


    # --- 4. Third Priority: IGR with Gene Name Extraction ---
    if "IGR:" in annotation_string:
        # Regex for IGR is different as the gene name follows 'before/after'
        match_igr_gene = re.search(r"(?:before|after)\s+([\w\d]+)(?:\(|\+|-)", annotation_string, re.IGNORECASE)
        
        if match_igr_gene:
            igr_gene_name = match_igr_gene.group(1)
            return f"IGR ({igr_gene_name})"
        else:
            return "IGR (Unknown Gene)"


    # --- 5. Intron Calculation (Fallback) ---
    if "Intron of" in annotation_string:
        # Note: We already extracted the gene name above, but the output here is the intron number.
        
        # Pattern for 'after exon X' (Intron X)
        match_after = re.search(r"after exon (\d+)", annotation_string, re.IGNORECASE)
        if match_after:
            #return f"Intron {int(match_after.group(1))} ({gene_name})"
            return f"Intron {int(match_after.group(1))}"

        # Pattern for 'before exon Y' (Intron Y-1)
        match_before = re.search(r"before exon (\d+)", annotation_string, re.IGNORECASE)
        if match_before:
            exon_number = int(match_before.group(1))
            if exon_number > 1:
                #return f"Intron {exon_number - 1} ({gene_name})"
                return f"Intron {exon_number - 1}"

    # --- 6. Default/Unclassified ---
    return annotation_string

def group_sv(df):
    """
    Groups rows based on ['Gene_St', 'Gene_End', 'Loc1'], treating each NaN/string 'nan' as unique.
    GroupID is assigned sequentially (0, 1, 2, ...) by factoring the combined group key 
    after the DataFrame is strictly sorted by ['#CHROM', 'POS'].
    """
    # 1. Define the grouping and sorting columns
    grouping_cols = ['Gene_St', 'Gene_End', 'Loc1','Loc2']
    
    # 2. Create a copy and prepare for strict VCF-style sorting
    df_copy = df.copy()

    # --- Robust VCF Sorting Logic using Categoricals ---
    
    # 2a. Convert POS to numeric
    df_copy.loc[:, 'POS'] = pd.to_numeric(df_copy['POS'], errors='coerce')
    
    # 2b. Create a numerically sortable list of chromosome names
    def numeric_chrom_key(chrom):
        chrom_str = str(chrom).lower().replace('chr', '')
        if chrom_str.isdigit():
            return int(chrom_str)
        elif chrom_str == 'x':
            return 23
        elif chrom_str == 'y':
            return 24
        elif chrom_str == 'm' or chrom_str == 'mt':
            return 25
        else:
            return 100 # Stable, high number for other names

    # Apply the key function and create a categorical sort key
    df_copy.loc[:, 'CHROM_KEY_NUM'] = df_copy['#CHROM'].apply(numeric_chrom_key)

    # 2c. Sort the DataFrame. This is CRUCIAL.
    strict_sort_cols = ['CHROM_KEY_NUM', 'POS']
    
    # The DataFrame is now strictly sorted by coordinate.
    df_copy = df_copy.sort_values(
        by=strict_sort_cols, 
        ascending=True, 
        na_position='last'
    ).reset_index(drop=True)
    
    df_copy = df_copy.drop(columns=['CHROM_KEY_NUM']) # Remove the temporary sort key

    # 3. Create unique, robust temporary grouping keys (to keep NaNs/string 'nan' unique)
    temp_grouping_cols = [f'temp_group_key_{col}' for col in grouping_cols]
    
    for col in grouping_cols:
        temp_col_name = f'temp_group_key_{col}'
        
        # Identify true NaNs and the string 'nan'
        is_missing = df_copy[col].isna() | (df_copy[col].astype(str).str.lower() == 'nan')
        
        # Convert all values to string (including valid ones)
        df_copy.loc[:, temp_col_name] = df_copy[col].astype(str)
        
        # Assign a unique string based on the sorted row index (0, 1, 2, 3...)
        if is_missing.any():
             unique_nan_ids = 'MISSING_UNIQUE_' + df_copy.index[is_missing].astype(str)
             df_copy.loc[is_missing, temp_col_name] = unique_nan_ids

    # 4. GUARANTEE Sequential Group ID Assignment using pd.factorize()
    # Combine the temporary grouping keys into a single string series
    combined_key = df_copy[temp_grouping_cols].astype(str).agg('_'.join, axis=1)

    # pd.factorize assigns numeric codes (starting at 0) in the order they appear.
    # Since the DataFrame is sorted, the codes will be assigned sequentially by coordinate.
    codes, unique_groups = pd.factorize(combined_key)
    df_copy.loc[:, 'GroupID'] = codes

    # 5. Calculate Count (using value_counts of the combined key)
    group_counts_df = combined_key.value_counts().reset_index()
    group_counts_df.columns = ['combined_key', 'Count']
    
    # Map the count back to the original DataFrame
    # Create the combined key on the main dataframe for merging
    df_copy.loc[:, 'combined_key'] = combined_key
    
    df_final = df_copy.merge(group_counts_df, on='combined_key', how='left')
    df_final = df_final.drop(columns=['combined_key'])

    # 6. Clean up temporary columns
    df_final = df_final.drop(columns=temp_grouping_cols)

    # 7. Reorder columns (The DataFrame remains sorted by #CHROM and POS)
    cols = [col for col in df_final.columns if col not in ['Count', 'GroupID']]
    df_final = df_final[cols + ['Count', 'GroupID']]
    df_final["BAR"] = df_final["BAR"].astype('float')
    df_final["UBAR"] = df_final["UBAR"].astype('int')
    df_final = df_final.sort_values(by=['GroupID', 'BAR', 'UBAR'], ascending=[True, False, False]).reset_index(drop=True)
    return df_final


def add_group_id(df1,df2):
	maxid=max(df1["GroupID"].to_list())+1
	df2["GroupID"]=df2.apply(lambda x:x["GroupID"]+maxid,axis=1)
	df=pd.concat([df1,df2],ignore_index=True)
	return df

def addvaf(df,BAM):
        AVG_INSERT_SIZE = 400
        STDEV_INSERT_SIZE = 2
        df["vaf"]=df.apply(lambda x: splitreads1.getVAF(BAM, x["#CHROM"],x["POS"],x["SV_type"],x["endchr"],x["endpos"]), axis=1)
        #df["vaf"]=df["vaf"].str.split("/").str[0]
        df["Total_Reads"]=df["vaf"].str.split("/").str[1]
        df["vaf"]=df["vaf"].str.split("/").str[0]
        return df

def svtype(x):
	if x["#CHROM"]!=x["endchr"]:
		return "TRA"
	if x["str1"]==x["str2"]:
		return "INV"
	else:
		return "INS"

def check_fusion(x):
	if x["site1"] is None or x["site2"] is None:
		return "-"
	if (x["#CHROM"]!=x["endchr"]) and ("IGR:" not in str(x["site1"]) and "IGR:" not in str(x["site2"])):
		return "Y"
	return "-"

df["Loc1"]=df.apply(lambda x: getloc(x["site1"]), axis=1)
df["Loc2"]=df.apply(lambda x: getloc(x["site2"]), axis=1)
df["chr1"]="chr"+df["chr1"].astype('str')
df["chr2"]="chr"+df["chr2"].astype('str')
df=df.rename(columns={"chr1":"#CHROM","chr2":"endchr","pos1":"POS","pos2":"endpos","gene1":"Gene_St","gene2":"Gene_End"})
#print (df)



extdf=pd.read_csv(sys.argv[2],sep="\t")
extdf["endchr"] = "chr"+extdf["ALT"].str.extract(r"chr(.*?):")
extdf["endpos"] = extdf["ALT"].str.extract(r":(\d+)").astype('int')
#print (extdf)

df=extdf.merge(df,how='left',on=["#CHROM","POS","endchr","endpos"])
df["APID"]=df.apply(lambda x: x["ID"].split("_")[0],axis=1)
df["BAR"]=df["FNUMS"].str.split(":").str[5]
df["UBAR"]=df["FNUMS"].str.split(":").str[6]
df=df.sort_values(by=["#CHROM","POS"])
print (df)
#df.to_csv("TB-fomatted.tsv",sep="\t",index=False)

df1=df.loc[df["ID"].str.contains("_1")]
df2=df.loc[df["ID"].str.contains("_2")]

df1= group_sv(df1)
df2= group_sv(df2)

df=add_group_id(df1,df2)
#df["Fusion"]=np.where( (df["Gene_St"].notna()) &  (df["Gene_End"].notna()) & (df["Gene_St"] != df["Gene_End"]) & (df["site1"]), "Y", "-")
df["GA_GB"]=df.apply(lambda x: check_fusion(x),axis=1)
df["SV_type"]=df.apply(lambda x: svtype(x), axis=1)
df["Repeat"] = (df["repName-repClass-repFamily:-site1"].fillna("-").astype(str) + "/" +df["repName-repClass-repFamily:-site2"].fillna("-").astype(str))
df["Loc1"]=(df["site1"].fillna("-").astype(str) + "/" +df["site2"].fillna("-").astype(str))
#df["SV_length"]=df.apply(lambda x: 0 if x["#CHROM"]!==x["endchr"] else abs(int(x["POS"])-int(x["endpos"])),axis=1)

df["SV_length"] = ((df["POS"].astype(int) - df["endpos"].astype(int)).abs() * (df["#CHROM"] == df["endchr"]))

df=addvaf(df,sys.argv[3])
df=df.rename(columns={"fusion":"Fusion"})
print (df)
print (list(df))
df.to_csv(sys.argv[4],sep="\t",index=False)

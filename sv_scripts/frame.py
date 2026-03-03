import pandas as pd
#from canonical_trans import *
import canonical_trans as ct

class CDSExon:
    """Represents a Coding Sequence (CDS) exon parsed from the GTF/DataFrame."""
    def __init__(self, start, end, phase, strand, exon_num, gene_id="N/A"):
        # Coordinates are 1-based, inclusive
        self.start = start 
        self.end = end
        self.phase = int(phase) # Column 'frame' in your DataFrame
        self.strand = strand
        self.gene_id = gene_id
        self.length = end - start + 1 # Genomic length
        self.exon_num = int(exon_num)
        
    def __repr__(self):
        return (f"CDSExon(start={self.start}, end={self.end}, "
                f"phase={self.phase}, exon_num={self.exon_num})")

def get_coords_for_intronic_break(breakpoint, exon_list_A, exon_list_B, strand_A, strand_B):
    # Gene A (Upstream): Find the last full exon BEFORE the break
    if strand_A == '+':
        # Filter for exons that end before the break, sort by end position descending
        valid_exons = [e for e in exon_list_A if e.end < breakpoint]
        last_exon = valid_exons[-1] # The last one biologically
        coord_A = last_exon.end # Snap to end of this exon
    else: # Strand -
        # Filter for exons that start after the break (genomically higher), sort ascending
        valid_exons = [e for e in exon_list_A if e.start > breakpoint]
        last_exon = valid_exons[0] # The 'first' genomically is the 'last' biologically
        coord_A = last_exon.start # Snap to start of this exon

    # Gene B (Downstream): Find the first full exon AFTER the break
    if strand_B == '+':
        # Filter exons starting after break
        valid_exons = [e for e in exon_list_B if e.start > breakpoint]
        first_exon = valid_exons[0]
        coord_B = first_exon.start # Snap to start
    else: # Strand -
        # Filter exons ending before break
        valid_exons = [e for e in exon_list_B if e.end < breakpoint]
        first_exon = valid_exons[-1]
        coord_B = first_exon.end # Snap to end
        
    return coord_A, coord_B

def read_gtf(path):
    df=pd.read_csv(path,sep="\t",names=["chr","source","feature","start","end","score","strand","frame","attribute"])
    df["gene_id"] = df["attribute"].str.extract(r'gene_id\s+"([^"]+)"', expand=False)
    df["transcript_id"] = df["attribute"].str.extract(r'transcript_id\s+"([^"]+)"', expand=False)
    df["gene_name"] = df["attribute"].str.extract(r'gene_name\s+"([^"]+)"', expand=False)
    df["exon_num"] = df["attribute"].str.extract(r'exon_number\s+"([^"]+)"', expand=False)
    #print (df)
    return df

def convert_dataframe_to_cds_exons(df):
    df = df.sort_values(by=['start', 'exon_num'], ascending=True).reset_index(drop=True)

    if df.empty:
        return [], None
    gene_strand = df['strand'].iloc[0]
 
    # 3. Apply conversion logic to each row
    cds_exons_list = []
    
    # Note: We use .itertuples() for generally faster row iteration than .iterrows()
    for row in df.itertuples():
        exon = CDSExon(
            start=row.start, 
            end=row.end, 
            phase=row.frame, 
            strand=row.strand, 
            exon_num=row.exon_num
        )
        cds_exons_list.append(exon)
        
    return cds_exons_list, gene_strand

def exons(hgdf,hgpath,gene,trans):
	df=hgdf.loc[(hgdf["transcript_id"]==trans) & (hgdf["gene_name"]==gene)]
	if df.empty:
		df=hgdf.loc[hgdf["transcript_id"]==trans]
		if df.empty:
			trans=ct.driver(hgpath,gene)
			print ("cannical", trans)
			if trans:
				df=hgdf.loc[hgdf["transcript_id"]==trans]
			else:
				return None,None
	df=df.loc[df["feature"]=="CDS"]
	if df.empty:
		return None,None
	df=df[["start","end","frame","strand","exon_num"]]
	#df=list(df.itertuples(index=False))
	cds_exons, strand = convert_dataframe_to_cds_exons(df)
	#print (cds_exons, strand)
	return cds_exons, strand

def check_frame_status(cds_length_A, phase_B):
    """
    Calculates if the fusion is in-frame based on Gene A's remainder 
    and Gene B's required starting phase.
    """
    # 1. Gene A contribution: How many bases remain from the last codon?
    remainder_A = cds_length_A % 3
    
    # 2. Gene B requirement: How many bases are needed to start a full codon?
    needed_by_B = (3 - phase_B) % 3
    
    if remainder_A == needed_by_B:
        return "IN-FRAME", 0
    else:
        # Calculate the magnitude of the frame shift
        shift = (remainder_A - needed_by_B) % 3
        return f"OUT-OF-FRAME (Shift of {shift} bp)", shift

def get_gene_A_contribution(breakpoint_A, cds_exons_A, strand_A):
    """Calculates the total coding length (L_A) retained from Gene A."""
    total_length = 0
    
    # Ensure exons are sorted by genomic coordinate (needed for logic to work)
    cds_exons_A.sort(key=lambda e: e.start)
    
    for exon in cds_exons_A:
        # Case 1: Exon is ENTIRELY KEPT (before the break)
        is_kept_full = (strand_A == '+' and exon.end < breakpoint_A) or \
                       (strand_A == '-' and exon.start > breakpoint_A)
        
        if is_kept_full:
            total_length += exon.length
            continue
        
        # Case 2: Break is INSIDE the coding exon (Partial Exon Kept)
        is_break_inside = (strand_A == '+' and exon.start <= breakpoint_A < exon.end) or \
                          (strand_A == '-' and exon.end >= breakpoint_A > exon.start)
        
        if is_break_inside:
            # Calculate partial length kept
            if strand_A == '+':
                # 5' partner on +, keeps the left side (up to breakpoint)
                partial_length = breakpoint_A - exon.start 
            else: # 5' partner on -, keeps the right side (from breakpoint)
                partial_length = exon.end - breakpoint_A
            
            total_length += partial_length
            
            # This is the last contribution from Gene A
            return total_length, exon.exon_num, "Partial"
            
        # Case 3: Exon is entirely past the break or in an intron after the break
        # In the intronic case, the loop terminates when the last full exon is passed.
        # If the break was intronic, total_length is the sum of full exons kept.
        # We return here as no more coding sequence is expected from A.
        if (strand_A == '+' and exon.start > breakpoint_A) or \
           (strand_A == '-' and exon.end < breakpoint_A):
            # Break must have been intronic, and we are on the first exon past the break
            return total_length, cds_exons_A[cds_exons_A.index(exon)-1].exon_num, "Full"

    # If the break was after the last exon, all CDS length is kept
    return total_length, cds_exons_A[-1].exon_num, "Full"

def get_gene_B_phase_and_exon(breakpoint_B, cds_exons_B, strand_B):
    """Determines the effective Phase (Frame) and the starting exon of Gene B."""
    
    # Ensure exons are sorted by genomic coordinate (needed for logic to work)
    cds_exons_B.sort(key=lambda e: e.start)
    
    for exon in cds_exons_B:
        # Case 1: Break is in an INT RON (We snap to the next full exon)
        # Gene B starts with a full, unsplit exon.
        is_intron_break = (strand_B == '+' and breakpoint_B < exon.start) or \
                          (strand_B == '-' and breakpoint_B > exon.end)
        
        if is_intron_break:
            return exon, exon.phase, "Full Exon Start"

        # Case 2: Break is in an EXON (Mid-CDS break)
        is_exon_break = (strand_B == '+' and exon.start <= breakpoint_B < exon.end) or \
                        (strand_B == '-' and exon.end >= breakpoint_B > exon.start)
        
        if is_exon_break:
            # Calculate bases skipped in this exon before the fusion starts
            if strand_B == '+':
                # 3' partner on +, starts from the base AFTER the breakpoint
                bases_skipped = breakpoint_B - exon.start + 1
            else: # 3' partner on -, starts from the base BEFORE the breakpoint
                bases_skipped = exon.end - breakpoint_B 
            
            # Effective Phase = (Original Phase + Skipped Bases) mod 3
            effective_phase = (exon.phase + bases_skipped) % 3
            return exon, effective_phase, "Partial Exon Start"

    return None, None, "No CDS Found"

def check_fusion_frame(exon_list_A, exon_list_B, strand_A, strand_B, pos1, pos2, geneA, geneB):
    """
    Checks the in-frame/out-of-frame status of a gene fusion.

    Args:
        exon_list_A/B (list): List of CDSExon objects for Gene A/B.
        strand_A/B (str): Strand ('+' or '-') of Gene A/B.
        pos1 (int): Genomic breakpoint position in Gene A.
        pos2 (int): Genomic breakpoint position in Gene B.
        geneA/B (str): Gene names.

    Returns:
        tuple: (Status, Details)
    """
    
    # --- Step 1: Calculate Gene A Contribution (L_A) ---
    L_A, exon_A_num, a_type = get_gene_A_contribution(pos1, exon_list_A, strand_A)
    
    if L_A == 0:
        return "NON-CODING (A starts after fusion point or is UTR only)", None

    # --- Step 2: Determine Gene B Start Phase (Phase_B) ---
    exon_B, phase_B, b_type = get_gene_B_phase_and_exon(pos2, exon_list_B, strand_B)
    
    if phase_B is None:
        return "NON-CODING (B starts after fusion point or is UTR/3' end only)", None

    # --- Step 3: Final Frame Check ---
    status, shift = check_frame_status(L_A, phase_B)
    
    # --- Step 4: Formatting Output ---
    details = (
        f"A ({geneA}): Contributes L_A={L_A} bp (Remainder: {L_A % 3}) up to Exon {exon_A_num} ({a_type}).\n"
        f"B ({geneB}): Starts at Exon {exon_B.exon_num} (Type: {b_type}, Phase: {phase_B}, Needs: {(3 - phase_B) % 3})."
    )
    
    return status, details

def driver_frame(hgdf,hgpath,geneA,transA,geneB,transB,pos1,pos2):
	print ("analyzing", geneA,transA,geneB,transB,pos1,pos2)
	exon_list_A, strand_A=exons(hgdf,hgpath,geneA,transA)
	exon_list_B, strand_B=exons(hgdf,hgpath,geneB,transB)
	if exon_list_A and exon_list_B:
		status, details = check_fusion_frame(exon_list_A, exon_list_B, strand_A, strand_B, pos1, pos2, geneA, geneB)
	else:
		status="No CDS"
		details=""
	print (status,details)
	return status



if __name__ == "__main__":
	df=pd.read_csv("orient.tsv",sep="\t")
	hgdf=read_gtf("hg19.refGene.gtf")
	df["status"]=df.apply(lambda x: driver_frame(hgdf,x["Gene_St"],x["transcript1"],x["Gene_End"],x["transcript2"],int(x["POS"]),int(x["endpos"])), axis=1)
	#driver(hgdf,geneA,transA,geneB,transB,pos1,pos2)
	print (df)
	df.to_csv("out.tsv",sep="\t")

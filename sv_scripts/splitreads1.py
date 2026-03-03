import pysam
import os # Added for path checking in a real-world scenario

def count_sv_metrics(bam_file, chrom, position, sv_type, partner_chrom=None):
    """
    Calculates the metrics needed for Structural Variant Allele Frequency (VAF), 
    customized for a specific SV type (TRA, INV, DEL, or INS).
    
    Args:
        bam_file (str): Path to the BAM file.
        chrom (str): Chromosome name (e.g., 'chr1').
        position (int): 1-based start position (the breakpoint) to check.
        sv_type (str): The type of structural variant ('TRA', 'INV', 'DEL', 'INS').
        partner_chrom (str, optional): The chromosome name of the translocation partner 
                                       (required for 'TRA'). Defaults to None.
        expected_insert_size (int): The mean insert size of the library.
        insert_size_stddev (int): The standard deviation of the insert size.
        stddev_factor (int): How many standard deviations to use for the cut-off.
            
    Returns:
        dict: A dictionary with counts.
    """
    expected_insert_size=500
    insert_size_stddev=50
    stddev_factor=3
    # Calculate insert size cut-offs (used for DEL/INS)
    #print (bam_file, chrom, position, sv_type, partner_chrom)
    #print (expected_insert_size,insert_size_stddev, stddev_factor)
    MAX_NORMAL_INSERT = expected_insert_size + (insert_size_stddev * stddev_factor)
    # Validate SV type and partner_chrom for Translocations
    sv_type = sv_type.upper() if isinstance(sv_type, str) else "TRA"
    if sv_type == 'TRA' and partner_chrom is None:
        return {'error': "A 'partner_chrom' must be provided for Translocation ('TRA')."}

    # Initialize counts
    discordant_count = 0
    split_count = 0
    total_coverage = 0
    
    # Set to track reads already counted as variant-supporting to avoid double-counting
    counted_sv_qnames = set() 
    
    # Open the BAM file
    try:
        if not os.path.exists(bam_file):
            return {'error': f"BAM file not found at {bam_file}. Check the path."}
            
        samfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        return {'error': f"Failed to open BAM file: {e}"}

    # Define a small, centered window around the breakpoint
    window_size = 80
    start = max(1, position - window_size)
    end = position + window_size
    counted_fragments_for_coverage = set()
    
    # --- Step 1: Count Variant Reads and Total Coverage ---
    # Fetch reads around the primary breakpoint
    for read in samfile.fetch(chrom, start, end):
        
        # Filter for primary, mapped reads
        if read.is_unmapped or read.is_secondary or read.is_duplicate:
            continue

        # Count towards total coverage (only one read per fragment/query_name)
        # Note: This logic assumes paired-end data. For single-end, simplify the logic.
        if read.is_paired and read.is_read1 and read.query_name not in counted_fragments_for_coverage:
             total_coverage += 1
             counted_fragments_for_coverage.add(read.query_name)
        elif not read.is_paired:
             total_coverage += 1

        # Skip read if it has already been counted as an SV-supporting read
        if read.query_name in counted_sv_qnames:
            continue

        # A. Split Read Check (Supplementary Alignment) - High-confidence SV support
        if read.has_tag('SA'):
            split_count += 1
            counted_sv_qnames.add(read.query_name)
            continue 

        # B. Discordant Read Pair Check (Read Pair Signal)
        # Check if paired and mapped, but NOT a 'proper pair'
        if read.is_paired and not read.mate_is_unmapped and read.is_proper_pair is False:
            
            # --- TRA/INV/DEL/INS-SPECIFIC DISCORDANT READ LOGIC ---
            is_sv_discordant = False
            
            if sv_type == 'TRA' and partner_chrom:
                # Translocation: Mate maps to a different, specified chromosome
                if read.next_reference_name == partner_chrom:
                    is_sv_discordant = True
                    
            elif sv_type == 'INV':
                # Inversion: Mate maps to the same chromosome but with incorrect orientation.
                # Expected orientation (standard FFPE): read.is_reverse != read.mate_is_reverse
                # Inversion signals: Both forward (-> ->) OR Both reverse (<- <-)
                if read.reference_name == read.next_reference_name:
                    if read.is_reverse == read.mate_is_reverse:
                        is_sv_discordant = True
         
            #elif sv_type == 'DEL':
                # Deletion: Insert size is much larger than expected
            #    if read.reference_name == read.next_reference_name:
            #        if abs(read.template_length) > MAX_NORMAL_INSERT:
            #            is_sv_discordant = True

            #elif sv_type == 'INS':
            else:
                # Insertion: Insert size is much smaller than expected (often near zero or negative TLEN)
                # However, many pipelines use large TLEN for INS as well. 
                # For simplicity here, we'll keep the large TLEN check from the original code
                # as a generic discordant check if it's the only one you're looking for.
                # A more precise INS check looks for very small/negative TLEN.
                if read.reference_name == read.next_reference_name:
                    if abs(read.template_length) > MAX_NORMAL_INSERT:
                        is_sv_discordant = True

            # Increment count if a specific SV-type discordant read is found
            if is_sv_discordant:
                discordant_count += 1
                counted_sv_qnames.add(read.query_name)
                    
    samfile.close()
    
    # --- Step 2: Calculate VAF ---
    total_sv_reads = discordant_count + split_count
    
    if total_coverage == 0:
        vaf = 0.0
    else:
        # VAF is the ratio of SV-supporting reads to total reads
        vaf = total_sv_reads / total_coverage
 
    resdict = {
        'discordant_reads': discordant_count,
        'split_reads': split_count,
        'total_sv_reads': total_sv_reads,
        'total_coverage': total_coverage,
        'VAF': round(vaf, 4),
        'sv_type_analyzed': sv_type
    }
    #print (resdict)
    return (resdict)
    #return vaf

def getVAF(BAM, CHROMOSOME, BREAKPOINT_POS, sv_type, partner_chrom, END_POS):
	bp1=count_sv_metrics(BAM, CHROMOSOME, int(BREAKPOINT_POS), sv_type, partner_chrom)
	bp2=count_sv_metrics(BAM, partner_chrom, int(END_POS), sv_type, CHROMOSOME)
	sv_reads=bp1["discordant_reads"]+bp2["discordant_reads"]+bp1["split_reads"]+bp2["discordant_reads"]
	tot_cov=bp1["total_coverage"]+bp2["total_coverage"]
	#print (bp1["discordant_reads"],bp2["discordant_reads"],bp1["split_reads"],bp2["discordant_reads"])
	#print (bp1["total_coverage"],bp2["total_coverage"])
	#print (bp1)
	#print (bp2)
	vaf=sv_reads/float(tot_cov) if sv_reads<tot_cov else 1
	#print (vaf)
	return str(vaf)+"/"+str(tot_cov)

# -------------------------------------------------------------
# CORRECTED EXAMPLE USAGE ⚠️
# -------------------------------------------------------------
# NOTE: The path '/Volumes/Molecular_Diagnostics/...' was updated for a Linux/cluster environment 
# but the local file system access will fail on any machine that does not have that path mounted.
# I'll use a placeholder for the final call.

if __name__ == '__main__':
	BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-22-17692_S4/TB-22-17692_S4_tumor.bam"
	CHROMOSOME = "chr11"
	BREAKPOINT_POS = 3748513
	sv_type="TRA"
	partner_chrom="chr2"
	END_POS=178165527

	# Call the correct function name with the defined variables
	#results = count_sv_metrics(BAM, CHROMOSOME, BREAKPOINT_POS, sv_type, partner_chrom, END_POS, AVG_INSERT_SIZE, STDEV_INSERT_SIZE)
	results= getVAF(BAM, CHROMOSOME, BREAKPOINT_POS, sv_type, partner_chrom, END_POS)
	print (results)

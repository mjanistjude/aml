import pysam
import os # Added for path checking in a real-world scenario

def count_sv_metrics(bam_file, chrom, position, sv_type, partner_chrom=None,
                      window_size=80, min_mapq=0):
    """
    Calculates the metrics needed for Structural Variant Allele Frequency (VAF), 
    customized for a specific SV type (TRA, INV, DEL, DUP, or INS).
    
    Args:
        bam_file (str): Path to the BAM file.
        chrom (str): Chromosome name (e.g., 'chr1').
        position (int): 1-based start position (the breakpoint) to check.
        sv_type (str): The type of structural variant ('TRA', 'INV', 'DEL', 'DUP', 'INS').
        partner_chrom (str, optional): The chromosome name of the translocation partner 
                                       (required for 'TRA'). Defaults to None.
        window_size (int): Half-width of the fetch window around `position`. Wider
                            windows are recommended for insert-size-based SV types
                            (DEL/INS) so that discordant pairs with large template
                            lengths still have their leftmost read fetched.
        min_mapq (int): Minimum mapping quality required to keep a read.
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

    # Define a centered window around the breakpoint (caller-configurable;
    # widen this for DEL/INS/TRA where discordant pairs can have large
    # template lengths and the leftmost read may sit outside a narrow window)
    start = max(1, position - window_size)
    end = position + window_size
    counted_fragments_for_coverage = set()
    
    # --- Step 1: Count Variant Reads and Total Coverage ---
    # Fetch reads around the primary breakpoint
    for read in samfile.fetch(chrom, start, end):
        
        # Filter for primary, mapped, QC-passing reads above the MAPQ threshold
        if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail:
            continue
        if read.mapping_quality < min_mapq:
            continue

        # Does this read's aligned (non-soft-clipped) span actually cross
        # the exact breakpoint base? This is what IGV's per-base depth
        # reflects.
        covers_breakpoint = (
            read.reference_end is not None
            and read.reference_start <= (position - 1) < read.reference_end
        )

        # --- Determine SV evidence status BEFORE deciding coverage ---
        # A split read is, by definition, soft-clipped right at the
        # junction - its aligned span often does NOT cross `position`
        # (covers_breakpoint can be False for genuine split evidence).
        # If coverage were gated on covers_breakpoint alone, such reads
        # would count toward split_reads/discordant_reads but never
        # toward total_coverage, reproducing the split_reads >
        # total_coverage problem from a different angle. So: figure out
        # whether this read qualifies as split/discordant evidence FIRST,
        # then let that also count it toward coverage below - a read
        # that IS the evidence is by construction an observation at this
        # locus, whether or not its CIGAR-mapped span reaches the exact
        # base.
        already_counted_as_sv = read.query_name in counted_sv_qnames
        is_split_evidence = (not already_counted_as_sv) and read.has_tag('SA')

        is_discordant_evidence = False
        if not already_counted_as_sv and not is_split_evidence:
            if read.is_paired and not read.mate_is_unmapped and read.is_proper_pair is False:

                # --- TRA/INV/DEL/DUP/INS-SPECIFIC DISCORDANT READ LOGIC ---
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

                elif sv_type == 'DEL':
                    # Deletion: Insert size is much larger than expected
                    if read.reference_name == read.next_reference_name:
                        if abs(read.template_length) > MAX_NORMAL_INSERT:
                            is_sv_discordant = True

                elif sv_type == 'DUP':
                    # Tandem duplication: read pair orientation is "outward
                    # facing" (RF) rather than the normal "inward facing" (FR)
                    # orientation - i.e. same opposite-strand pairing as a
                    # normal pair, but with the forward/reverse mates' relative
                    # order flipped: the forward-strand mate maps DOWNSTREAM of
                    # its reverse-strand mate (or equivalently, the
                    # reverse-strand mate maps UPSTREAM of its forward mate).
                    # This is evaluated purely from each read's own alignment
                    # flags/position - not from the caller's SV-level STRANDS/
                    # str1,str2 annotation, which only classifies the overall
                    # call and says nothing about any individual read's
                    # orientation.
                    if read.reference_name == read.next_reference_name:
                        if read.is_reverse != read.mate_is_reverse:
                            if (not read.is_reverse and read.reference_start > read.next_reference_start) or \
                               (read.is_reverse and read.reference_start < read.next_reference_start):
                                is_sv_discordant = True

                elif sv_type == 'INS':
                    # Insertion: for a pure read-pair signal, effective insert size
                    # is often smaller than expected (or unremarkable) since the
                    # inserted sequence isn't in the reference; INS support mainly
                    # comes from split/soft-clipped reads (handled above via 'SA').
                    # This branch is intentionally conservative and will likely
                    # undercount INS-supporting discordant pairs - split reads are
                    # the primary signal for insertions.
                    if read.reference_name == read.next_reference_name:
                        if abs(read.template_length) < (expected_insert_size - insert_size_stddev * stddev_factor):
                            is_sv_discordant = True

                else:
                    # Fallback for any other/unrecognized sv_type (e.g. generic
                    # 'FUSION' label): treat abnormally large insert size as the
                    # discordant signal. TRA, INV, DEL, DUP and INS all now have
                    # their own branches above - this only catches values that
                    # don't match any of those exactly (e.g. typos, or a label
                    # from an unexpected caller).
                    if read.reference_name == read.next_reference_name:
                        if abs(read.template_length) > MAX_NORMAL_INSERT:
                            is_sv_discordant = True

                is_discordant_evidence = is_sv_discordant

        # Count towards total coverage: one count per unique READ (not
        # per fragment) that EITHER spans the exact breakpoint base OR is
        # itself being counted as split/discordant evidence for this
        # breakpoint (see comment above).
        #
        # FIXED: deduping by query_name alone collapsed R1 and R2 of the
        # same fragment into a single coverage count. IGV's per-base
        # depth (like samtools depth/pileup) counts each aligned READ
        # separately - if both mates independently span the breakpoint
        # (common with short templates, e.g. cfDNA, where mates overlap),
        # that fragment should contribute 2 to coverage, not 1. Keying on
        # (query_name, is_read1) restores that: R1 and R2 are counted
        # separately, while a single mate's own primary + supplementary
        # (chimeric) alignment records - which share the same is_read1 -
        # still collapse to one count, so a split read's two segments
        # aren't double-counted as two units of depth for one physical
        # read.
        coverage_key = (read.query_name, read.is_read1 if read.is_paired else None)
        if (covers_breakpoint or is_split_evidence or is_discordant_evidence) \
           and coverage_key not in counted_fragments_for_coverage:
            total_coverage += 1
            counted_fragments_for_coverage.add(coverage_key)

        # Commit the SV-evidence counts determined above.
        if is_split_evidence:
            split_count += 1
            counted_sv_qnames.add(read.query_name)
        elif is_discordant_evidence:
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

def getVAF(BAM, CHROMOSOME, BREAKPOINT_POS, sv_type, partner_chrom, END_POS,
           window_size=80, min_mapq=0):
	bp1=count_sv_metrics(BAM, CHROMOSOME, int(BREAKPOINT_POS), sv_type, partner_chrom,
	                     window_size=window_size, min_mapq=min_mapq)
	bp2=count_sv_metrics(BAM, partner_chrom, int(END_POS), sv_type, CHROMOSOME,
	                     window_size=window_size, min_mapq=min_mapq)

	# Surface breakpoint-level errors (e.g. bad BAM path) instead of silently
	# proceeding with a KeyError or nonsense numbers.
	if 'error' in bp1:
		return {'error': f"breakpoint 1: {bp1['error']}"}
	if 'error' in bp2:
		return {'error': f"breakpoint 2: {bp2['error']}"}

	# FIXED: previously double-counted bp2's discordant reads and never
	# included bp2's split reads at all.
	sv_reads = bp1["discordant_reads"] + bp2["discordant_reads"] + bp1["split_reads"] + bp2["split_reads"]
	tot_cov = bp1["total_coverage"] + bp2["total_coverage"]

	cov_bp1 = bp1["total_coverage"]
	cov_bp2 = bp2["total_coverage"]

	sv_reads_bp1 = bp1["discordant_reads"] + bp1["split_reads"]
	sv_reads_bp2 = bp2["discordant_reads"] + bp2["split_reads"]

	def _safe_vaf(reads, cov):
		# Same clamping logic as the combined VAF below, applied per-breakpoint.
		if cov == 0:
			return 0.0
		if reads >= cov:
			return 1.0
		return round(reads / float(cov), 4)

	vaf_bp1 = _safe_vaf(sv_reads_bp1, cov_bp1)
	vaf_bp2 = _safe_vaf(sv_reads_bp2, cov_bp2)

	print ("BP1 Dis:",bp1["discordant_reads"],"BP2 Dis:",bp2["discordant_reads"],"BP1 Split:",bp1["split_reads"],"BP2 Split:",bp2["split_reads"],"total_split_reads:",int(bp1["split_reads"])+int(bp2["split_reads"]))
	print ("BP1 tot:",cov_bp1,"BP2 tot:",cov_bp2)
	#print (bp1)
	#print (bp2)

	if tot_cov == 0:
		return {'error': 'no coverage at either breakpoint'}

	if sv_reads >= tot_cov:
		# This should not normally happen; flag it instead of silently
		# clamping to 1.0, since it usually indicates overlapping windows
		# or a double-counting bug rather than a true 100% VAF.
		print(f"WARNING: sv_reads ({sv_reads}) >= tot_cov ({tot_cov}); "
		      f"clamping combined VAF to 1.0. Check for window overlap or double-counting.")
		vaf = 1.0
	else:
		vaf = round(sv_reads / float(tot_cov), 4)

	ratioA = bp1["split_reads"] / cov_bp1 if cov_bp1 else 0
	ratioB = bp2["split_reads"] / cov_bp2 if cov_bp2 else 0

	return {
		'VAF': vaf,
		'total_coverage': tot_cov,
		'coverage_bp1': cov_bp1,
		'coverage_bp2': cov_bp2,
		'VAF_bp1': vaf_bp1,
		'VAF_bp2': vaf_bp2,
		'split_reads_total': bp1["split_reads"] + bp2["split_reads"],
		'split_reads_bp1': bp1["split_reads"],
		'split_reads_bp2': bp2["split_reads"],
		'discordant_reads_bp1': bp1["discordant_reads"],
		'discordant_reads_bp2': bp2["discordant_reads"],
		'ratioA':ratioA,
		'ratioB':ratioB,
	}

# -------------------------------------------------------------
# CORRECTED EXAMPLE USAGE ⚠️
# -------------------------------------------------------------
# NOTE: The path '/Volumes/Molecular_Diagnostics/...' was updated for a Linux/cluster environment 
# but the local file system access will fail on any machine that does not have that path mounted.
# I'll use a placeholder for the final call.

if __name__ == '__main__':
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1/combined/combined/Output/varcall/19-364-0009-0007/19-364-0009-0007_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1_Dx/combined/Output/SD-19-00575_S3/SD-19-00575_S3_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1_Dx/combined/Output/SD-19-00787_S4/SD-19-00787_S4_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1_Dx/combined/Output/23P-073MP0024_S2/23P-073MP0024_S2_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1_Dx/combined/Output/23P-151MP0036_S6/23P-151MP0036_S6_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1_Dx/combined/Output/23P-073MP0024_S2/23P-073MP0024_S2_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260626_LH00399_0044_B237JYLLT3/combined/Output/varcall/24P-144BR0043-144IP0009/24P-144BR0043-144IP0009_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1/combined/combined/Output/varcall/19-364-0009-0007/19-364-0009-0007_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260626_LH00399_0044_B237JYLLT3/combined/Output/varcall/24P-079MP0013/24P-079MP0013_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260626_LH00399_0044_B237JYLLT3/combined/Output/varcall/24P-137IP0010/24P-137IP0010_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260626_LH00399_0044_B237JYLLT3/combined/Output/varcall/24P-158MP0019/24P-158MP0019_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260626_LH00399_0044_B237JYLLT3/combined/Output/varcall/24P-205MP0031/24P-205MP0031_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/additional_samples/20260324_LH00399_0040_A22FLC2LT1/diagnosis_samples/Output/24P-130MP08_S1/24P-130MP08_S1_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/additional_samples/20260324_LH00399_0040_A22FLC2LT1/diagnosis_samples/Output/25P-MP6044_S3/25P-MP6044_S3_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/additional_samples/20260324_LH00399_0040_A22FLC2LT1/diagnosis_samples/Output/26P-MP0354_S2/26P-MP0354_S2_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML032027_D1_S99_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-20-07289/TB-20-07289_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML032043_D1_S66_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML031434_D1_S144_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-19-14332_S11/TB-19-14332_S11_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-21-44403_S12/TB-21-44403_S12_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-21-44395/TB-21-44395_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML032598_D1_S117_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML032732_D1_S55_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-21-47489/TB-21-47489_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-21-47420_S7/TB-21-47420_S7_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML03
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-22-19582/TB-22-19582_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-22-19392_S10/TB-22-19392_S10_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML033661_D1_S135_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-22-16473/TB-22-16473_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-22-18538/TB-22-18538_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260306_LH00399_0036_A23LKMGLT4/Output/varcall/TB-22-16415_S5/TB-22-16415_S5_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260306_LH00399_0036_A23LKMGLT4/Output/varcall/TB-22-18637_S1/TB-22-18637_S1_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML033284_D1_S125_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-22-12218/TB-22-12218_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260306_LH00399_0036_A23LKMGLT4/Output/varcall/TB-22-12192_S2/TB-22-12192_S2_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML033421_D1_S128_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-22-11805/TB-22-11805_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-22-13970/TB-22-13970_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-22-14324_S1/TB-22-14324_S1_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML031473_D1_S35_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-19-13804/TB-19-13804_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-19-15771/TB-19-15771_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-19-13805_S6/TB-19-13805_S6_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-19-15753_S2/TB-19-15753_S2_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML031855_D1_S102_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML033573_D1_S133_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML032879_D1_S121_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260310_LH00399_0037_B235M3MLT3/combined/Output/varcall/TB-21-49220/TB-21-49220_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-21-51508/TB-21-51508_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML031206_D1_S90_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML032193_D1_S54_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-20-09211_S14/TB-20-09211_S14_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-20-18749_S11/TB-20-18749_S11_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-20-09218/TB-20-09218_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-20-18844/TB-20-18844_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/cohort/diagnosis_samples/Output/SJAML032773_D1_S63_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-21-45253/TB-21-45253_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-21-47477/TB-21-47477_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260306_LH00399_0036_A23LKMGLT4/Output/varcall/TB-21-45433_S8/TB-21-45433_S8_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260306_LH00399_0036_A23LKMGLT4/Output/varcall/TB-21-47457_S4/TB-21-47457_S4_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1/combined/combined/Output/varcall/23P-109IP0011/23P-109IP0011_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1/combined/combined/Output/varcall/20-037-0632/20-037-0632_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1/combined/combined/Output/varcall/20-055-1160/20-055-1160_tumor.bam"
	BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1/combined/combined/Output/varcall/24P-082MP0018/24P-082MP0018_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1/combined/combined/Output/varcall/24P-103IP0003-103BR0008/24P-103IP0003-103BR0008_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/additional_samples/genomic_samples/Output/varcall/26P-MP0001218/26P-MP0001218_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/additional_samples/genomic_samples/Output/varcall/25P-MP0006372/25P-MP0006372_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260626_LH00399_0044_B237JYLLT3/combined/Output/varcall/26P-MP0002481/26P-MP0002481_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260626_LH00399_0044_B237JYLLT3/combined/Output/varcall/26P-MP0001098/26P-MP0001098_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-18-11511/TB-18-11511_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-18-12935/TB-18-12935_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-20-05562/TB-20-05562_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-20-05635_S3/TB-20-05635_S3_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-20-07259_S9/TB-20-07259_S9_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-20-06211/TB-20-06211_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-20-07503/TB-20-07503_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-20-06187_S2/TB-20-06187_S2_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-20-07695_S8/TB-20-07695_S8_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-19-10302_S5/TB-19-10302_S5_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-19-14056/TB-19-14056_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-19-10008/TB-19-10008_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-21-42750/TB-21-42750_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-21-42930_S6/TB-21-42930_S6_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-21-45626_S1/TB-21-45626_S1_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-21-45657/TB-21-45657_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20251010_LH00399_0029_A237NN7LT4/combined/Output/varcall/TB-22-17693/TB-22-17693_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/Output/varcall/TB-22-17692_S4/TB-22-17692_S4_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-22-09925/TB-22-09925_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260306_LH00399_0036_A23LKMGLT4/Output/varcall/TB-22-09918_S6/TB-22-09918_S6_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-22-11601_S5/TB-22-11601_S5_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-20-02247_S7/TB-20-02247_S7_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-20-02140/TB-20-02140_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-20-03923/TB-20-03923_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260310_LH00399_0037_B235M3MLT3/combined/Output/varcall/TB-22-16990/TB-22-16990_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260306_LH00399_0036_A23LKMGLT4/combined/Output/varcall/TB-22-15042/TB-22-15042_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-22-15282_S8/TB-22-15282_S8_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-22-16968_S4/TB-22-16968_S4_tumor.bam"
	#BAM="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/20260310_LH00399_0037_B235M3MLT3/Output/varcall/TB-21-51481_S9/TB-21-51481_S9_tumor.bam"

	#chr11:118353021::chr10:21972689
	#chr10:22248346::chr11:118357562
	sv_type="TRA"  # FIXED: was "FUSION", which silently skipped the
	               # different-chromosome discordant-read logic and only
	               # counted split reads for this translocation

	#basepath="/clinical/ccs01/dept/PATH/Molecular_Diagnostics/Eval/MRD_AML16_workinggroup/data_analysis/patient_cohort/genomic_dna/20260701_LH00399_0047_B22FTVTLT1_Dx/combined/Output/"
	#samples=["23P-109IP0011","24P-010IP0021","24P-086IP0009","19-364-0009-0007","20-037-0632","20-055-1160","24P-082MP0018","24P-103IP0003-103BR0008","24P-144BR0043-144IP0009"]
	#samples=["SD-19-00575"]
	#strs=["chr11:118354770::chr6:168251544","chr11:118354770::chr6:168251544","chr11:118354770::chr6:168251544","chr11:118357568::chr10:21958529","chr11:118357568::chr10:21958529","chr11:118357568::chr10:21958529","chr11:118355803::chr19::18618042","chr11:118355803::chr19::18618042","chr11:118355803::chr19::18618042"]
	#codict={}
	#for s,c in zip(samples,strs):
	#	codict[s]=c

	svstr="chr11:118354770::chr6:168251544"
	svstr="chr11:118357568::chr10:21958529"
	svstr="chr11:118355803::chr19:18618042"
	svstr="chr11:118354705::chr9:20358799"
	svstr="chr12:12019386::chr13:28608347"
	svstr="chr11:118359879::chr19:18562125"
	svstr="chr11:3761995::chr7:27208111"
	svstr="chr11:118352529::chr11:118342963"
	svstr="chr11:118354579::chrX:118826006"
	svstr="chr5:176660736::chr11:3759537"
	svstr="chr5:176660736::chr11:3759537"
	svstr="chr21:36225137::chr8:93082219"
	svstr="chr11:118353152::chr19:6277614"
	svstr="chr11:85673816::chr10:21870890"
	svstr="chr21:36221466::chr8:93061387"
	svstr="chr21:36228985::chr8:93080602"
	svstr="chr11:118359356::chr10:21910454"
	svstr="chr11:118353278::chr10:21958239"
	svstr="chr16:67131617::chr16:15815286"
	svstr="chr16:15815295::chr16:67120179"
	svstr="chr11:3755036::chr5:133863705"
	svstr="chr1:3076800::chr11:122931393"
	svstr="chr12:11832670::chr7:156763115"
	svstr="chr12:11832670::chr7:156763115"
	svstr="chr11:118355621::chr9:20376086"
	svstr="chr11:118353021::chr10:21972689"
	svstr="chr11:118352529::chr11:118342963"
	svstr="chr11:118354579::chrX:118826006"
	svstr="chr11:118352529::chr11:118342963"
	svstr="chr13:28608256::chr13:28608300"
	svstr="chr11:3755036::chr5:133863705"
	svstr="chr11:85673816::chr10:21870890"
	svstr="chr11:118355803::chr19:18618042"
	if 1:
		CHROMOSOME=svstr.split("::")[0].split(":")[0]
		BREAKPOINT_POS = int(svstr.split("::")[0].split(":")[1])
		partner_chrom = svstr.split("::")[1].split(":")[0]
		END_POS= int(svstr.split("::")[1].split(":")[1])

		# Call the correct function name with the defined variables
		#results = count_sv_metrics(BAM, CHROMOSOME, BREAKPOINT_POS, sv_type, partner_chrom, END_POS, AVG_INSERT_SIZE, STDEV_INSERT_SIZE)
		results= getVAF(BAM, CHROMOSOME, BREAKPOINT_POS, sv_type, partner_chrom, END_POS)
		print (results)

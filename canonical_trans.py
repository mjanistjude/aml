import sys
import csv
from collections import defaultdict, namedtuple

GTFRecord = namedtuple("GTFRecord", ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attrs"])
ASSUME_EXONS_ARE_CDS_FOR_NM = True   # If CDS features are missing, treat NM_ exons as CDS (approximate)
PERMISSIVE_TRANSCRIPT_FALLBACK = True  # If preferred transcript not found, choose canonical by gene

def parse_gtf_attributes(attr_str):
    """
    Parse GTF attributes: key "value"; -> dict.
    """
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, val = part.split(" ", 1)
            attrs[key] = val.strip().strip('"')
    return attrs

def read_gtf(gtf_path):
    """
    Build transcript models from GTF: exons and CDS (if present).
    Returns dict: transcripts[transcript_id] = {
        'chr': str,
        'strand': '+'|'-',
        'gene_name': str or None,
        'exons': [(start, end), ...] in transcription order,
        'cds':   [(start, end), ...] in transcription order,
    }
    """
    transcripts = {}
    exons = defaultdict(list)
    cds = defaultdict(list)
    t_strand = {}
    t_chr = {}
    t_gene = {}

    with open(gtf_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            rec = GTFRecord(
                chr=fields[0],
                source=fields[1],
                feature=fields[2],
                start=int(fields[3]),
                end=int(fields[4]),
                score=fields[5],
                strand=fields[6],
                frame=fields[7],
                attrs=parse_gtf_attributes(fields[8]),
            )
            tid = rec.attrs.get("transcript_id") or rec.attrs.get("transcriptId") or rec.attrs.get("transcript")
            if not tid:
                continue

            gene_name = rec.attrs.get("gene_name") or rec.attrs.get("gene_id")
            if gene_name:
                t_gene[tid] = gene_name
            t_strand[tid] = rec.strand
            t_chr[tid] = rec.chr

            if rec.feature.lower() == "exon":
                exons[tid].append((rec.start, rec.end))
            elif rec.feature.lower() == "cds" or rec.feature.upper() == "CDS":
                cds[tid].append((rec.start, rec.end))

    for tid in set(list(exons.keys()) + list(cds.keys()) + list(t_strand.keys())):
        strand = t_strand.get(tid, "+")
        chrom = t_chr.get(tid, None)
        gene_name = t_gene.get(tid, None)
        exon_list = exons.get(tid, [])
        cds_list = cds.get(tid, [])

        # Sort intervals in transcription order
        if strand == "+":
            exon_list.sort(key=lambda x: (x[0], x[1]))
            cds_list.sort(key=lambda x: (x[0], x[1]))
        else:
            exon_list.sort(key=lambda x: (x[0], x[1]), reverse=True)
            cds_list.sort(key=lambda x: (x[0], x[1]), reverse=True)

        transcripts[tid] = {
            "chr": chrom,
            "strand": strand,
            "gene_name": gene_name,
            "exons": exon_list,
            "cds": cds_list,
        }
    return transcripts

def pick_canonical_transcript_by_gene(transcripts, gene_name):
    """
    Given a gene symbol, pick a canonical transcript:
      - Prefer NM_ transcripts (protein-coding).
      - Among candidates, choose longest CDS length; if no CDS, choose longest exon span.
    Returns (chosen_tid, chosen_transcript) or (None, None) if none found.
    """
    # Collect candidates
    candidates = []
    for tid, t in transcripts.items():
        if t.get("gene_name") == gene_name:
            cds_len = sum(e - s + 1 for s, e in t["cds"]) if t["cds"] else 0
            exon_span = sum(e - s + 1 for s, e in t["exons"]) if t["exons"] else 0
            candidates.append((tid, t, cds_len, exon_span))
    if not candidates:
        return None, None

    # Prefer NM_ candidates
    nm_candidates = [c for c in candidates if c[0].startswith("NM_")]
    pool = nm_candidates if nm_candidates else candidates

    # Choose longest CDS; if CDS length ties/zero, fall back to longest exon span
    best = max(pool, key=lambda x: (x[2], x[3]))
    return best[0], best[1]

def driver(path,gene):
	#transcripts=read_gtf("hg19.refGene.gtf")
	transcripts=read_gtf(path)
	tid, trans=pick_canonical_transcript_by_gene(transcripts,gene)
	#print (tid)
	#print (trans)
	return tid

if __name__ == "__main__":
	tid=driver("hg19.refGene.gtf","KMT2A")
	print (tid)

import sys
import re
from collections import namedtuple, defaultdict
import pysam
from multiprocessing import Pool
import random
import os
import subprocess


#read alignment constants
MIN_ALIGNED_LENGTH = 7000
MIN_ALIGNED_RATE = 0.9
MAX_SEGMENTS = 10
MIN_SEGMENT_MAPQ = 10
MIN_SEGMENT_LENGTH = 100


ReadSegment = namedtuple("ReadSegment", ["read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                                         "strand", "read_length", "haplotype", "mapq", "genome_id"])


#TODO: merge with get_segment below
ReadSegmentLegacy = namedtuple("ReadSegmentLegacy", ["read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                                                     "strand", "read_length", "mismatch_rate"])
def _get_segment_legacy(read_id, ref_id, ref_start, strand, cigar, num_mismatch):
    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    #read_end = 0

    for token in re.findall("[\d]{0,}[A-Z]{1}", cigar):
        op = token[-1]
        op_len = int(token[:-1])

        if op in "HS":
            if first_clip:
                read_start = op_len
            read_length += op_len
        first_clip = False

        if op in "M=X":
            read_aligned += op_len
            ref_aligned += op_len
            read_length += op_len
        if op == "D":
            ref_aligned += op_len
        if op == "I":
            read_aligned += op_len
            read_length += op_len

    ref_end = ref_start + ref_aligned
    read_end = read_start + read_aligned

    #mm_rate = num_mismatch / (read_aligned + 1)
    length_diff = abs(ref_aligned - read_aligned)
    mm_rate = (num_mismatch - length_diff) / (read_aligned + 1)

    #print(read_id, mm_rate, mm_rate_2)

    if strand == "-":
        read_start, read_end = read_length - read_end, read_length - read_start

    return ReadSegmentLegacy(read_start, read_end, ref_start, ref_end, read_id,
                       ref_id, strand, read_length, mm_rate)
###


def check_read_mapping_confidence(sam_text_entry, min_aln_length, min_aligned_rate,
                                  max_read_error, max_segments):
    fields = sam_text_entry.split()

    read_id, flags, chr_id, position = fields[0:4]
    cigar = fields[5]
    ref_id, ref_start = fields[2], int(fields[3])

    is_supplementary = int(flags) & 0x800
    is_secondary = int(flags) & 0x100
    is_unmapped = int(flags) & 0x4
    strand = "-" if int(flags) & 0x10 else "+"

    if is_unmapped:
        return False

    sa_tag = ""
    num_mismatches = 0
    #de_tag = None
    for tag in fields[11:]:
        if tag.startswith("SA"):
            sa_tag = tag[5:]
        if tag.startswith("NM"):
            num_mismatches = int(tag[5:])
        #if tag.startswith("de"):
        #    de_tag = float(tag[5:])

    segments = [_get_segment_legacy(read_id, ref_id, ref_start, strand, cigar, num_mismatches)]
    #print(de_tag, segments[0].mismatch_rate)
    if sa_tag:
        for sa_aln in sa_tag.split(";"):
            if sa_aln:
                sa_fields = sa_aln.split(",")
                sa_ref, sa_ref_pos, sa_strand, sa_cigar, sa_mismatches = \
                        sa_fields[0], int(sa_fields[1]), sa_fields[2], sa_fields[3], int(sa_fields[5])
                segments.append(_get_segment_legacy(read_id, sa_ref, sa_ref_pos, sa_strand, sa_cigar, sa_mismatches))
    segments.sort(key=lambda s: s.read_start)

    read_length = segments[0].read_length

    WND_LEN = 100
    read_coverage = [0 for x in range(read_length // WND_LEN)]
    weighted_mm_sum = 0
    total_segment_length = 0
    for seg in segments:
        for i in range(seg.read_start // WND_LEN, seg.read_end // WND_LEN):
            read_coverage[i] = 1
        weighted_mm_sum += seg.mismatch_rate * (seg.read_end - seg.read_start)
        total_segment_length += seg.read_end - seg.read_start
    aligned_length = sum(read_coverage) * WND_LEN
    mean_mm = weighted_mm_sum / (total_segment_length + 1)

    if (is_secondary or 
            is_unmapped or 
            aligned_length < min_aln_length or 
            aligned_length / read_length < min_aligned_rate or
            len(segments) > max_segments or
            mean_mm > max_read_error):
        return False

    else:
        return True
        #if len(segments) > 1:
        #    print(read_length, len(segments), aligned_length, aligned_length / read_length, mean_mm)
        #if is_supplementary:
        #    new_read_id = read_id + "_suppl_" + str(unique_id)
        #    unique_id += 1
        #    new_flags = int(flags) & ~0x800
        #    aln.query_name = new_read_id
        #    aln.flag = new_flags



cigar_parser = re.compile("[0-9]+[MIDNSHP=X]")
def get_segment(read_id, ref_id, ref_start, strand, cigar, haplotype, mapq, genome_id):
    """
    Parses cigar and generate ReadSegment structure with alignment coordinates
    """
    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    #read_end = 0

    for token in cigar_parser.findall(cigar):
        op = token[-1]
        op_len = int(token[:-1])

        if op == "H" or op == "S":
            if first_clip:
                read_start = op_len
            read_length += op_len
        first_clip = False

        if op == "M" or op == "=" or op == "X":
            read_aligned += op_len
            ref_aligned += op_len
            read_length += op_len
        if op == "D":
            ref_aligned += op_len
        if op == "I":
            read_aligned += op_len
            read_length += op_len

    ref_end = ref_start + ref_aligned
    read_end = read_start + read_aligned

    if strand == "-":
        read_start, read_end = read_length - read_end, read_length - read_start

    return ReadSegment(read_start, read_end, ref_start, ref_end, read_id,
                       ref_id, strand, read_length, haplotype, mapq, genome_id)


def get_split_reads(bam_file, ref_id, max_read_error, genome_id):
    """
    Yields set of split reads for each contig separately. Only reads primary alignments
    and infers the split reads from SA alignment tag
    """
    filtered_reads = set()
    alignments = []

    aln_file = pysam.AlignmentFile(bam_file, "rb")
    for aln in aln_file.fetch(ref_id, multiple_iterators=True):
        fields = aln.to_string().split()
        read_id, flags, chr_id, position = fields[0:4]
        cigar = fields[5]
        mapq = int(fields[4])
        ref_id, ref_start = fields[2], int(fields[3])

        is_supplementary = int(flags) & 0x800
        is_secondary = int(flags) & 0x100
        is_unmapped = int(flags) & 0x4
        is_primary = not (is_supplementary or is_secondary or is_unmapped)
        strand = "-" if int(flags) & 0x10 else "+"

        if not check_read_mapping_confidence(aln.to_string(), MIN_ALIGNED_LENGTH, MIN_ALIGNED_RATE,
                                             max_read_error, MAX_SEGMENTS):
            if is_primary:
                filtered_reads.add(aln.query_name)
            continue

        #if is_supplementary or is_secondary or is_unmapped:
        if is_secondary or is_unmapped:
            continue

        sa_tag = None
        hp_tag = None
        for tag in fields[11:]:
            if tag.startswith("SA"):
                sa_tag = tag[5:]
            if tag.startswith("HP"):
                hp_tag = tag[5:]
        if not hp_tag:
            haplotype = 0
        else:
            haplotype = int(hp_tag)

        new_segment = get_segment(read_id, ref_id, ref_start, strand, cigar, haplotype, mapq, genome_id)
        if new_segment.mapq >= MIN_SEGMENT_MAPQ and new_segment.read_end - new_segment.read_start >= MIN_SEGMENT_LENGTH:
            alignments.append(new_segment)

    return filtered_reads, alignments


def get_all_reads_parallel(bam_file, num_threads, aln_dump_stream, max_read_error, genome_id):

    all_reference_ids = [r for r in pysam.AlignmentFile(bam_file, "rb").references]
    random.shuffle(all_reference_ids)
    tasks = [(bam_file, r, max_read_error, genome_id) for r in all_reference_ids]

    parsing_results = None
    with Pool(num_threads) as p:
        parsing_results = p.starmap(get_split_reads, tasks)

    all_filtered_reads = set()
    segments_by_read = defaultdict(list)
    for filtered, alignments in parsing_results:
        all_filtered_reads |= filtered
        for aln in alignments:
            segments_by_read[aln.read_id].append(aln)

    all_reads = []
    for read in segments_by_read:
        if read not in all_filtered_reads:
            aln_dump_stream.write(str(read) + "\n")
            for seg in segments_by_read[read]:
                aln_dump_stream.write(str(seg) + "\n")
            aln_dump_stream.write("\n")

            segments = segments_by_read[read]
            segments.sort(key=lambda s: s.read_start)
            all_reads.append(segments)
        #else:
        #    four.write("Filtered: " + str(read))

    return all_reads


def get_segments_coverage(bam_file, segments, num_threads):
    print("Computing segments coverage", file=sys.stderr)

    random.shuffle(segments)
    tasks = [(bam_file, s[0], s[1], s[2]) for s in segments]

    seg_coverage = None
    with Pool(num_threads) as p:
        seg_coverage = p.starmap(_get_median_depth, tasks)

    coverage_map = {}
    for seg, coverage in zip(segments, seg_coverage):
        coverage_map[seg] = coverage

    return coverage_map


def _get_median_depth(bam_path, ref_id, ref_start, ref_end):
    SAMTOOLS_BIN = "samtools"
    samtools_out = subprocess.Popen("{0} coverage {1} -r '{2}:{3}-{4}' -q 10 -l 100"
                                     .format(SAMTOOLS_BIN, bam_path, ref_id, ref_start, ref_end),
                                    shell=True, stdout=subprocess.PIPE).stdout

    for line in samtools_out:
        if line.startswith(b"#"):
            continue
        fields = line.split()
        coverage = float(fields[6])
        return coverage

    return None

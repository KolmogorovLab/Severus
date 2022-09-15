#!/usr/bin/env python3

"""
Script to filter out possibly mismapped alignments in sam
"""

import os
import sys
import pysam
import re
import random
from collections import namedtuple, defaultdict
import subprocess
import signal
import multiprocessing as mp


ReadSegment = namedtuple("ReadSegment", ["read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                                         "strand", "read_length", "mismatch_rate"])


def get_segment(read_id, ref_id, ref_start, strand, cigar, num_mismatch):
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

    return ReadSegment(read_start, read_end, ref_start, ref_end, read_id,
                       ref_id, strand, read_length, mm_rate)


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

    segments = [get_segment(read_id, ref_id, ref_start, strand, cigar, num_mismatches)]
    #print(de_tag, segments[0].mismatch_rate)
    if sa_tag:
        for sa_aln in sa_tag.split(";"):
            if sa_aln:
                sa_fields = sa_aln.split(",")
                sa_ref, sa_ref_pos, sa_strand, sa_cigar, sa_mismatches = \
                        sa_fields[0], int(sa_fields[1]), sa_fields[2], sa_fields[3], int(sa_fields[5])
                segments.append(get_segment(read_id, sa_ref, sa_ref_pos, sa_strand, sa_cigar, sa_mismatches))
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


def filter_alignments(bam_in, bam_out, contig_ids, min_aligned_length, max_read_error, bam_index):
    #MIN_ALIGNED_LENGTH = 10000
    MAX_SEGMENTS = 3
    MIN_ALIGNED_RATE = 0.9
    #MAX_READ_ERROR = 0.1

    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index)
    bam_writer = pysam.AlignmentFile(bam_out, "wb", template=bam_reader)
    for ctg in contig_ids:
        for aln in bam_reader.fetch(ctg):
            line = aln.to_string()
            if check_read_mapping_confidence(line, min_aligned_length, MIN_ALIGNED_RATE, max_read_error, MAX_SEGMENTS):
                bam_writer.write(aln)


def filter_alignments_parallel(bam_in, bam_out, num_threads, min_aligned_length, max_read_error,
                               bam_index=None):
    all_reference_ids = [r for r in pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).references]
    random.shuffle(all_reference_ids)

    bams_to_merge = []
    threads = []
    chunk_size = len(all_reference_ids) // num_threads + 1

    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    for i in range(num_threads):
        contigs_list = all_reference_ids[i * chunk_size : (i + 1) * chunk_size]
        if not contigs_list:
            continue

        bam_out_part = bam_out + "_part_" + str(i)
        bams_to_merge.append(bam_out_part)
        threads.append(mp.Process(target=filter_alignments, args=(bam_in, bam_out_part, contigs_list,
                                                                  min_aligned_length, max_read_error,
                                                                  bam_index)))

    signal.signal(signal.SIGINT, orig_sigint)

    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
            if t.exitcode != 0:
                raise Exception("One of the processes exited with code: {0}".format(t.exitcode))
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()
        raise

    pysam.merge("-@", str(num_threads), bam_out, *bams_to_merge)
    for bam in bams_to_merge:
        os.remove(bam)


def main():
    filter_alignments_parallel(sys.argv[1], sys.argv[2], 16, min_aligned_length=10000, max_read_error=0.1)


if __name__ == "__main__":
    main()



#!/usr/bin/env python

from __future__ import division
import sys
import re


def fasta_lengths(handle):
    records = {}
    record_row = re.compile("^>")
    current_id = None
    current_seq = 0
    for line in handle:
        if record_row.match(line) is not None:
            if current_id is not None:
                records[current_id] = current_seq
            id_ = line.strip(">").split(" ")[0].rstrip()
            current_id = id_
            current_seq = 0
        else:
            current_seq += len(line.strip())

    return records


def chunk_seqs(lengths, n=10, penalty=15):

    lengths = sorted(lengths.items(), key=lambda t: t[1], reverse=True)

    total_length = sum(l[1] for l in lengths)
    av_chunk_length = total_length / n

    chunks = [[] for i in range(n)]

    for i in range(len(lengths)):
        sid, length = lengths[i]
        dists = []
        for j in range(n):
            chunk_length = sum(t[1] for t in chunks[j]) + length
            diff = av_chunk_length - chunk_length
            # If projected length is greater than average, penalise it
            if diff < 0:
                diff *= -1
                diff *= penalty
            dists.append(diff)

        dist_min = min(dists)
        dist_min_idx = [j for j, d in enumerate(dists) if d == dist_min][0]
        chunks[dist_min_idx].append(lengths[i])

    return chunks

lengths = fasta_lengths(sys.stdin)

try:
    n = int(sys.argv[1])
except:
    n = 12

chunks = chunk_seqs(lengths, n=n, penalty=15)

for chunk in chunks:
    seqids = [t[0] for t in chunk]
    print(" ".join(seqids))

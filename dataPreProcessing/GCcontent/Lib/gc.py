"""
Functions to compute GC and CpG content of a DNA sequence.
"""

import sys
import pandas as pd
import numpy as np
from numpy.testing import assert_array_equal
import gzip

def read_fasta(name):
    """
    Read an .fa.gz file in as an np.array.
    This function assumes a single chromosome per file.
    """
    if name.endswith('gz'):
        with gzip.open(name, "rb") as f:
            lines = f.read().split("\n")
    else:
        with open(name, "r") as f:
            lines = f.read().split("\n")
    header = lines[0].replace(">", "")
    seq = np.array([letter.lower() for letter in list("".join(lines[1:]))])
    assert seq.dtype == np.dtype("S1") # Make sure the datatype is a single character.
    return seq

def compute_gc(seq, window = 101):
    "Compute gc content for the given DNA sequence seq."
    gc = np.logical_or(seq == "g", seq == "c").astype(int) # 1's for gc, 0's for ta
    ser = pd.Series(gc)
    ma = ser.rolling(window = window, center = True, min_periods = 0).mean()
    return ma.values

def compute_gc_slow(seq, window = 101):
    """
    Slow way to compute gc content. For unit tests, to make sure the fast way is
    correct.
    """
    def gc_subseq(subseq):
        return np.mean(np.logical_or(subseq == "g", subseq == "c"))
    n = len(seq)
    width = (window - 1) / 2
    gcs = np.empty(n)
    def in_start_fringe(i): return i < width
    def in_end_fringe(i): return i > n - width
    for i in range(n):
        start = max(0, i-width)
        end = min(n, i + width + 1) # Add 1 since slicing exclusive at right endpoint.
        subseq = seq[start:end]
        gcs[i] = gc_subseq(subseq)
    return gcs

def compute_cpg(seq, window = 201):
    "Compute CpG content for the given DNA sequence"
    ind = cpg_indicator(seq)
    ser = pd.Series(ind)
    ma = ser.rolling(window = window, center = True, min_periods = 0).mean()
    return ma.values

def cpg_indicator(seq):
    """
    Replace the sequence of nucleotides by a sequence that is 1 for CG pairs,
    else 0.
    For example: attacggacgg becomes 00001100110
    """
    cpg_start = np.logical_and(seq[:-1] == "c", seq[1:] == "g").astype(int)
    return np.logical_or(np.concatenate([np.array([0]), cpg_start]),
                         np.concatenate([cpg_start, np.array([0])]))

################################################################################


if __name__ == "__main__":
    infile = sys.argv[1] # fasta file
    outfile = sys.argv[2] 
    window_size = int(sys.argv[3])
    seq = read_fasta(infile)
    gc = compute_gc(seq, window=window_size)
    np.savetxt(outfile,gc,fmt="%.4f",delimiter=' ')



################################################################################

# Unit tests.

def test_gc():
    "Check the gc content computation by 2 different methods."
    fname = "/net/noble/vol1/data/ucsc/goldenPath/hg19/chromosomes/chr10.fa.gz"
    seq = read_fasta(fname)
    subseq = np.concatenate([seq[:100], seq[100000:101000], seq[-100:]])
    gc_fast = compute_gc(subseq)
    gc_slow = compute_gc_slow(subseq)
    assert_array_equal(gc_fast, gc_slow)

def test_cpg_indicator():
    "Make sure that the cpg indicator function does the right thing."
    seq = np.array(list("attacggacgg"))
    expected = np.array([0,0,0,0,1,1,0,0,1,1,0])
    assert_array_equal(expected, cpg_indicator(seq))

# Alon Diament, Tuller Lab, June 2022.

from typing import Iterable
import re
import numpy as np


the_code = {
    'AAA': 'K',
    'AAC': 'N',
    'AAG': 'K',
    'AAT': 'N',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'AGA': 'R',
    'AGC': 'S',
    'AGG': 'R',
    'AGT': 'S',
    'ATA': 'I',
    'ATC': 'I',
    'ATG': 'M',
    'ATT': 'I',
    'CAA': 'Q',
    'CAC': 'H',
    'CAG': 'Q',
    'CAT': 'H',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'GAA': 'E',
    'GAC': 'D',
    'GAG': 'E',
    'GAT': 'D',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'TAA': '*',
    'TAC': 'Y',
    'TAG': '*',
    'TAT': 'Y',
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TGA': '*',
    'TGC': 'C',
    'TGG': 'W',
    'TGT': 'C',
    'TTA': 'L',
    'TTC': 'F',
    'TTG': 'L',
    'TTT': 'F',
}
nt2codon_dict = {k: chr(i+1) for i, k in enumerate(the_code.keys())}
valid_re = re.compile('[^ACGTacgt]+')


def nt2aa(seq_nt, validate_seq=True):
    if is_str_iter(seq_nt):
        # returning an empty string for bad seqs to preserve indexing
        seq_aa = [nt2aa(s)
                  if not validate_seq or
                  is_valid_seq(s) else ''
                  for s in seq_nt]
        print(f'validate_seq: ignored {sum([len(s) == 0 for s in seq_aa])} ambiguous sequences.')
        return seq_aa

    n = len(seq_nt)
    n = n - n % 3  # ignore partial codons
    seq_nt = seq_nt.upper().replace('U', 'T')
    seq_aa = ''.join([the_code[seq_nt[i:i+3]]
                      if seq_nt[i:i+3] in the_code
                      else 'X'
                      for i in range(0, n, 3)])

    n_nt = len(seq_nt) / 3
    n_aa = len(seq_aa)
    if n_aa < n_nt:
        print(f'nt2aa: ignored {n_nt-n_aa} ambiguous AA.')

    return seq_aa


def nt2codon(seq_nt, validate_seq=True):
    """ converts a sequence in NT alphabet and returns a sequence in codon
        alphabet, which is defined as follows:
        the index of the character at a position equals the index of the codon in
        the sorted list of all 64 triplets.
        this allows for all string algorithms, including the Chimera approach, to
        work seemlessly.
        ignores partial codons.
        Alon Diament, Tuller Lab, August 2018 (MATLAB), June 2022 (Python). """

    if is_str_iter(seq_nt):
        # returning an empty string for bad seqs to preserve indexing
        seq_cod = [nt2codon(s)
                   if not validate_seq or
                   is_valid_seq(s) else ''
                   for s in seq_nt]
        print(f'validate_seq: ignored {sum([len(s) == 0 for s in seq_cod])} ambiguous sequences.')
        return seq_cod

    if not len(seq_nt):
        return ''

    n = len(seq_nt)
    n = n - n % 3  # ignore partial codons
    seq_nt = seq_nt.upper().replace('U', 'T')
    seq_cod = ''.join([nt2codon_dict[seq_nt[i:i+3]]
                       if seq_nt[i:i+3] in nt2codon_dict
                       else chr(0)
                       for i in range(0, n, 3)])

    n_nt = len(seq_nt) / 3
    n_cod = len(seq_cod)
    if n_cod < n_nt:
        print(f'nt2codon: ignored {n_nt-n_cod} ambiguous codons.')

    return seq_cod


def codon2nt(seq_cod):
    if is_str_iter(seq_cod):
        return [codon2nt(s) for s in seq_cod]

    if not len(seq_cod):
        return ''

    codon_list = ['NNN'] + list(nt2codon_dict.keys())

    return ''.join([codon_list[ord(c)] for c in seq_cod])


def sample_seqs_from_blocks(all_blocks, n, random_seed=42):
    """
    generate multiple synonymous sequences by sampling different
    combinations of blocks from `all_blocks`. return `n` sampled
    sequences.
    """
    np.random.seed(random_seed)

    seqs = []
    for _ in range(n):
        seqs.append(''.join([np.random.choice(b, 1)[0]
                             for b in all_blocks]))

    return seqs


def rand_seq(n):
    return ''.join(np.random.choice(['A', 'C', 'G', 'T'], size=n))


def compare_seq(seq1, seq2):
    """ list edits between 2 identically long sequences.
    """
    return [[i, s1, s2] for i, (s1, s2) in
            enumerate(zip(seq1, seq2)) if s1 != s2]


def is_str_iter(obj):
    return not isinstance(obj, str) and isinstance(obj, Iterable)


def is_valid_seq(seq):
    return valid_re.search(seq) is None

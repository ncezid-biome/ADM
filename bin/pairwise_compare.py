#!/usr/bin/env python3
'''
import os
import argparse
import glob
from pathlib import Path
from Bio import SeqIO
import utilities

def check_diff_by_primer(seq1, seq2):
    if seq1 == seq2 or seq1 == utilities.revcomp(seq2):
        return False
    return True

# def compare(query_idx, target_idx, primer_list, numeric_flag, diff_only):
#     diff_count = 0
#     total_primer = len(primer_list)
    
#     for primer in primer_list:
#         row_seq_list = [val for key, val in query_idx.items() if primer in key]
#         col_seq_list = [val for key, val in target_idx.items() if primer in key]

#         if len(row_seq_list) != 1 or len(col_seq_list) != 1:
#             total_primer -= 1
#             continue

#         if check_diff_by_primer(row_seq_list[0].seq, col_seq_list[0].seq):
#             diff_count += 1

#     if numeric_flag:
#         return diff_count / total_primer if total_primer else "NA"
#     elif diff_only:
#         return diff_count
#     else:
#         return f"{diff_count}/{total_primer}"

def compare(query_idx, target_idx, primer_list, numeric_flag, diff_only):
    diff_count = 0
    total_primer = 0

    for primer in primer_list:
        query_recs = query_idx.get(primer, [])
        target_recs = target_idx.get(primer, [])

        if len(query_recs) != 1 or len(target_recs) != 1:
            continue  # skip if not exactly 1 vs 1 amplimer match
        total_primer += 1

        if check_diff_by_primer(query_recs[0].seq, target_recs[0].seq):
            diff_count += 1

    if total_primer == 0:
        return "NA"

    if numeric_flag:
        return diff_count / total_primer
    elif diff_only:
        return diff_count
    else:
        return f"{diff_count}/{total_primer}"


def build_primer_index(fasta_file, primer_list):
    primer_idx = {}

    for record in SeqIO.parse(fasta_file, 'fasta'):
        for primer in primer_list:
            if primer in record.id:
                if primer not in primer_idx:
                    primer_idx[primer] = []
                primer_idx[primer].append(record)
                # No break! because one record might match multiple primers (rare, but possible)

    return primer_idx



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', required=True, help='Query fasta file')
    # parser.add_argument('--all', nargs='+', required=True, help='All fasta files')
    parser.add_argument('--all', required=True, help='Folder for All fasta files')
    parser.add_argument('--primers', required=True, help='Oligos file')
    parser.add_argument('--output', required=True, help='Output row file')
    parser.add_argument('--numeric', action='store_true')
    parser.add_argument('--diff_only', action='store_true')
    args = parser.parse_args()

    primers = utilities.Primers(args.primers)
    primer_list = primers.pnames

    query_name = Path(args.query).stem
    # query_idx = SeqIO.to_dict(SeqIO.parse(args.query, 'fasta'))
    query_idx = build_primer_index(args.query, primer_list)


    # print (f'query_name is: {query_name}')
    # print (f'primer_list size: {len(primer_list)}')

    row = [query_name]
    # for target_fasta in sorted(args.all):
    # for target_fasta in sorted(glob.glob(f'{args.all}/**/*.fasta', recursive=True)):
    #     # target_idx = SeqIO.to_dict(SeqIO.parse(target_fasta, 'fasta'))
    #     target_idx = build_primer_index(target_fasta, primer_list)
    #     row.append(str(compare(query_idx, target_idx, primer_list, args.numeric, args.diff_only)))

    all_fastas = sorted(glob.glob(f'{args.all}/**/*.fasta', recursive=True))
    total = len(all_fastas)
    print(f"Comparing against {total} target FASTA files...")

    for i, target_fasta in enumerate(all_fastas):
        target_idx = build_primer_index(target_fasta, primer_list)
        result = compare(query_idx, target_idx, primer_list, args.numeric, args.diff_only)
        row.append(str(result))
        
        if (i + 1) % 100 == 0 or i == total - 1:
            print(f"Processed {i + 1}/{total} files")


    with open(args.output, 'w') as f:
        f.write(','.join(row) + '\n')

if __name__ == '__main__':
    main()


'''

import argparse
import glob
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import utilities
import hashlib
import os

def hash_sequence(sequence: str) -> str:

    md5 = hashlib.md5(sequence.encode("utf-8"))
    max_bits_in_result = 56
    p = (1 << max_bits_in_result) - 1
    rest = int(md5.hexdigest(), 16)
    result = 0
    while rest != 0:
        result = result ^ (rest & p)
        rest = rest >> max_bits_in_result
    return str(result)

hash_cache = {}
def get_hash(seq):
    seq_str = str(seq)
    if seq_str not in hash_cache:
        # hash_cache[seq_str] = hashlib.md5(seq_str.encode()).hexdigest()
        hash_cache[seq_str] = hash_sequence(seq_str)
    return hash_cache[seq_str]

# revcomp_cache = {}
# def revcomp_cached(seq):
#     if seq not in revcomp_cache:
#         revcomp_cache[seq] = utilities.revcomp(seq)
#     return revcomp_cache[seq]

# --- Cached reverse complement, now returns lower-case ---
revcomp_cache = {}
def revcomp_cached(seq):
    seq_upper = seq.upper()
    if seq_upper not in revcomp_cache:
        revcomp_cache[seq_upper] = utilities.revcomp(seq_upper)
    return revcomp_cache[seq_upper]

# def check_diff_by_primer(seq1, seq2):
#     if len(seq1) > 40:
#         h1 = get_hash(seq1)
#         h2 = get_hash(seq2)
#         # h2_rc = get_hash(utilities.revcomp(seq2))
#         h2_rc = get_hash(revcomp_cached(seq2))
#         return not (h1 == h2 or h1 == h2_rc)
#     else:
#         return not (seq1 == seq2 or seq1 == utilities.revcomp(seq2))

# --- Sequence comparison, now case-insensitive ---
def check_diff_by_primer(seq1, seq2):
    # convert sequences to lowercase
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()

    if len(seq1) > 40:
        h1 = get_hash(seq1)
        h2 = get_hash(seq2)
        h2_rc = get_hash(revcomp_cached(seq2))
        return not (h1 == h2 or h1 == h2_rc)
    else:
        return not (seq1 == seq2 or seq1 == revcomp_cached(seq2))
    
    
def compare(query_idx, target_idx, primer_list, numeric_flag, diff_only):
    diff_count = 0
    total_primer = 0

    for primer in primer_list:
        query_recs = query_idx.get(primer, [])
        target_recs = target_idx.get(primer, [])

        if len(query_recs) != 1 or len(target_recs) != 1:
            continue  # skip if not exactly 1 vs 1 amplimer match

        total_primer += 1
        if check_diff_by_primer(query_recs[0].seq, target_recs[0].seq):
            diff_count += 1

    if total_primer == 0:
        return "NA"

    if numeric_flag:
        return diff_count / total_primer
    elif diff_only:
        return diff_count
    else:
        return f"{diff_count}/{total_primer}"


def build_primer_index(fasta_file, primer_list):
    primer_idx = defaultdict(list)

    for record in SeqIO.parse(fasta_file, 'fasta'):
        for primer in primer_list:
            if primer in record.id:
                primer_idx[primer].append(record)
                break  # stop checking other primers for this record

    return primer_idx


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', required=True, help='Query fasta file')
    parser.add_argument('--all', required=True, help='Folder or single fasta file')
    parser.add_argument('--primers', required=True, help='Oligos file')
    parser.add_argument('--output', required=True, help='Output row file')
    parser.add_argument('--numeric', action='store_true')
    parser.add_argument('--diff_only', action='store_true')
    args = parser.parse_args()

    primers = utilities.Primers(args.primers)
    primer_list = primers.pnames

    query_name = Path(args.query).stem
    query_idx = build_primer_index(args.query, primer_list)

    # Load all target indexes just once
    # all_fastas = sorted(glob.glob(f'{args.all}/**/*.fasta', recursive=True))
    # all_indexes = [(Path(f).stem, build_primer_index(f, primer_list)) for f in all_fastas]
    
    # Determine if --all is folder or file
    if os.path.isdir(args.all):
        all_fastas = sorted(glob.glob(f'{args.all}/**/*.fasta', recursive=True))
    else:
        # If single file, use its parent directory
        all_fastas = sorted(glob.glob(f'{Path(args.all).parent}/**/*.fasta', recursive=True))

    all_indexes = [(Path(f).stem, build_primer_index(f, primer_list)) for f in all_fastas]

    row = [query_name]

    # for sample_name, target_idx in tqdm(all_indexes, desc="Comparing"):
    for sample_name, target_idx in all_indexes:
        result = compare(query_idx, target_idx, primer_list, args.numeric, args.diff_only)
        row.append(str(result))

    with open(args.output, 'w') as f:
        f.write(','.join(row) + '\n')


if __name__ == '__main__':
    main()

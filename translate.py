#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    rna_sequence = rna_sequence.upper()
    amino_acid_list = []
    while True:
        if len(rna_sequence) < 3:
            break
        codon = rna_sequence[0:3]
        remaining_seq = rna_sequence[3:]
        rna_sequence = remaining_seq
        aa = genetic_code[codon]
        if aa == "*":
            break
        amino_acid_list.append(aa)
    return "".join(amino_acid_list)

def get_all_translations(rna_sequence, genetic_code):
    rna_sequence = rna_sequence.upper()
    number_of_bases = len(rna_sequence)
    last_codon_index = number_of_bases - 3
    if last_codon_index < 0:
        return ('')
    amino_acid_seq_list = []
    for base_index in range(last_codon_index +1):
        codon = rna_sequence[base_index: base_index +3]
        if codon == "AUG":
            aa_seq = translate_sequence(
            rna_sequence = rna_sequence[base_index:],
            genetic_code = genetic_code)
            if aa_seq:
                amino_acid_seq_list.append(aa_seq)
    return amino_acid_seq_list

def get_reverse(sequence):
    if sequence:
        seq=sequence.upper()
        reverse_seq=seq[::-1]
        return reverse_seq
    else:
        return ''

def get_complement(sequence):
    if sequence:
        seq=list(sequence.upper())
        complements = {'C':'G','G':'C','U':'A','A':'U'}
        seq=[complements[base] for base in seq]
        return ''.join(seq)
    else:
        return ''    


def reverse_and_complement(sequence):
    if sequence:
        seq = get_reverse(sequence)
        rev_comp = get_complement(seq)
        return rev_comp
    else:
        return ''

def get_longest_peptide(rna_sequence, genetic_code):
    peptides = get_all_translations(rna_sequence,
        genetic_code = genetic_code)
    rev_comp_seq = reverse_and_complement(rna_sequence)
    rev_comp_peptides = get_all_translations(rna_sequence = rev_comp_seq,
        genetic_code = genetic_code)
    peptides += rev_comp_peptides
    if not peptides:
        return ""
    if len(peptides) <2:
        return peptides [0]
    most_number_of_bases = -1
    longest_peptide_index = -1
    for peptide_index, aa_seq in enumerate(peptides):
        if len(aa_seq) > most_number_of_bases:
            longest_peptide_index = peptide_index
            most_number_of_bases = len(aa_seq)
    return peptides[longest_peptide_index]

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")

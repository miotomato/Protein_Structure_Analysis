# Python 3.7

"""
File:           parsing.py
Author:         Shengyuan Wang
Date:           Jun 1, 2019
Last update:    May 26, 2020

Purpose:        Collect and parse individual protein individual chains data, including protein description,
                protein sequence with start/end position, associate nucleotide sequence with accession number and
                start/end position, DSSP secondary structure.

Input files:    PISCES id30 list, PISCES individual protein chain description and sequence, PDB files, DSSP files.
Output files:   PISCES clean file, TBLASTN files, clean TBLASTN files, final protein individual chain files with
                codon, amino acid, residue number, coordinates, DSSP.

Notices:        1. Proteins data are downloaded from PDB have the following limitations: Resolution <= 2.0, R-factor <=
                0.25, sequence length between 100-1000, sequence identity <= 30%.
                2. PDB protein sequence files need to be cleaned before use.
                3. After obtain TBLASTN files, check organism, remove '-', check 3' -> 5' or 5' -> 3',
                perform longest_common_sequence between PISCES protein sequence and TBLASTN nucleotide translated
                sequence.
                4. The longest_common_sequences in this case set to at least 30 residues.
"""

import os
import glob
import pandas as pd
from time import sleep

import shutil
import urllib.request as request
from contextlib import closing

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def prepare_pdb_file(input_dir, pdb_ori_id, pdb_ori_csv, output_dir, pdb_clean_file):
    """
    Protein chains are downloaded from PDB website with a customized report: PDB ID, Chain ID, R Free,
    Structure Title, Exp. Method, Resolution, Classification, Uniprot Acc, Source, EC No, Chain Length, Sequence and
    Secondary Structure in csv file.

    pdb_ori_id file contains protein id and chain number, pdb_ori_csv contains protein entries details. Both files
    are downloaded from PDB website with limits as follow: Resolution <= 2.0, R-factor free <= 0.25, seq_len between
    100 - 1000, chain sequence identity <= 30%.
    """
    f_id = list(open(os.path.join(input_dir, pdb_ori_id), 'r'))
    f_csv = pd.read_csv(os.path.join(input_dir, pdb_ori_csv))
    g = open(os.path.join(output_dir, pdb_clean_file), 'a')

    # Convert number to alphabet for protein chains.
    num_to_alphabet = {'1': 'A', '2': 'B', '3': 'C', '4': 'D', '5': 'E', '6': 'F', '7': 'G', '8': 'H', '9': 'I',
                       '10': 'J', '11': 'K', '12': 'L', '13': 'M', '14': 'N', '15': 'O', '16': 'P', '17': 'Q',
                       '18': 'R', '19': 'S', '20': 'T', '21': 'U', '22': 'V', '23': 'W', '24': 'X', '25': 'Y',
                       '26': 'Z'}

    pdbid_list = []
    for i in range(len(f_id)):
        pdbid_list.append(f_id[i][:4] + num_to_alphabet[f_id[i][5]])

    f_csv_id = f_csv['PDB ID']
    f_csv_chain = f_csv['Chain ID']
    f_csv_rf = f_csv['R Free']
    f_csv_title = f_csv['Structure Title']
    f_csv_r = f_csv['Resolution']
    f_csv_class = f_csv['Classification']
    f_csv_uni = f_csv['Uniprot Acc']
    f_csv_source = f_csv['Source']
    f_csv_len = f_csv['Chain Length']
    f_csv_ss = f_csv['Sequence and Secondary Structure']

    for i in range(len(f_csv_id)):
        if len(f_csv_id[i]) > 4:
            pid = f_csv_id[i][0] + f_csv_id[i][4] + f_csv_id[i][-2:]
        else:
            pid = f_csv_id[i]
        pdbid = str(pid + f_csv_chain[i])
        if pdbid in pdbid_list:
            loc = f_csv_ss[i].find('#')
            seq = str(f_csv_ss[i][:loc])
            ss = str(f_csv_ss[i][loc + 1:])
            g.writelines('>' + pdbid + '\t' + str("%.3f" % float(f_csv_rf[i])) + '\t' + str(f_csv_title[i]) + '\t' +
                         str(f_csv_r[i]) + '\t' + str(f_csv_class[i]) + '\t' + str(f_csv_uni[i]) + '\t' +
                         str(f_csv_source[i]) + '\t' + str(f_csv_len[i]) + '\n')
            g.writelines(seq + '\n')
            g.writelines(ss + '\n' + '\n')


def tblastn(input_dir, pdb_clean_file, output_dir):
    """
    Input protein file with description on the first line and FASTA sequence on the next line. Very time-consuming,
    better to perform multi-thread.
    """
    f = list(open(os.path.join(input_dir, pdb_clean_file), 'r'))

    for i in range(len(f)):
        if f[i][0] == '>':
            file = str(f[i] + f[i + 1])
            result_handle = NCBIWWW.qblast("tblastn", "nt", file)

            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            g = open(os.path.join(output_dir, "%s_tblastn" % str(f[i][1:6])), 'a')
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 0.1:
                        g.writelines('>' + alignment.title + '\n')
                        g.writelines(str(hsp.align_length) + '\t' + str(hsp.identities) + '\t' + str(hsp.query_start) +
                                     '\t' + str(hsp.query_end) + '\t' + str(hsp.query) + '\t' + str(hsp.sbjct_start) +
                                     '\t' + str(hsp.sbjct_end) + '\t' + str(hsp.sbjct) + '\t' + str(hsp.match) + '\n')
            g.close()


def find(str, ch):
    for i, ltr in enumerate(str):
        if ltr == ch:
            yield i


def extract_correct_tblastn(input_dir, pdb_clean_file, tblastn_dir, output):
    f = list(open(os.path.join(input_dir, pdb_clean_file), 'r'))
    g = open(os.path.join(input_dir, output), 'a')

    for i in range(len(f)):
        if f[i][0] == '>':
            pdb_id = f[i][1:6]
            pdb_organism = [m.split('\t', 7)[6] for m in [f[i][:-1]]][0]
            pdb_organism_short = pdb_organism.lower()
            ss = f[i + 2][:-1]

            if ' ' in pdb_organism:
                pdb_organism_split = pdb_organism.split()
                pdb_organism_short = pdb_organism_split[0].lower()

            # Human immunodeficiency virus was written short in TBLASTN files.
            if "Human immunodeficiency virus" in pdb_organism_short:
                pdb_organism_short = 'hiv'

            # Check if TBLASTN file exists in database.
            if os.path.isfile(os.path.join(tblastn_dir, "%s_tblastn" % pdb_id)):
                f_tblastn = list(open(os.path.join(tblastn_dir, "%s_tblastn" % pdb_id), 'r'))

                for j in range(len(f_tblastn)):
                    if f_tblastn[j][0] == '>':
                        try:
                            end = int(list(find(f_tblastn[j], " "))[2])
                        except IndexError:
                            end = len(f_tblastn[j])

                        # Remove predicted and synthetic proteins.
                        if 'predicted' in f_tblastn[j][:end].lower() or 'synthetic' in f_tblastn[j][:end].lower():
                            pass
                        else:
                            acc_start = int(list(find(f_tblastn[j], "|"))[2]) + 1
                            acc_end = int(list(find(f_tblastn[j], "|"))[3])
                            acc = f_tblastn[j][acc_start:acc_end]
                            align_len = [m.split('\t', 1)[0] for m in [f_tblastn[j + 1][:-1]]][0]
                            identity = [m.split('\t', 2)[1] for m in [f_tblastn[j + 1][:-1]]][0]
                            query_start = [m.split('\t', 3)[2] for m in [f_tblastn[j + 1][:-1]]][0]
                            query_end = [m.split('\t', 4)[3] for m in [f_tblastn[j + 1][:-1]]][0]
                            query = [m.split('\t', 5)[4] for m in [f_tblastn[j + 1][:-1]]][0]
                            subject_1 = int([m.split('\t', 6)[5] for m in [f_tblastn[j + 1][:-1]]][0])
                            subject_2 = int([m.split('\t', 7)[6] for m in [f_tblastn[j + 1][:-1]]][0])
                            subject_start, subject_end = str(min(subject_1, subject_2)), str(max(subject_1, subject_2))
                            subject = [m.split('\t', 8)[7] for m in [f_tblastn[j + 1][:-1]]][0]

                            query_ss = ss[int(query_start) - 1:int(query_end)]

                            # Only extract TBLASTN results if organism same as PDB or (match identity >=0.96 and
                            # match length >= 100).
                            if pdb_organism_short in f_tblastn[j][:end].lower() or \
                                    (int(align_len) >= 100 and int(identity) / int(align_len) >= 0.96):
                                g.writelines('>' + pdb_id + '\t' + acc + '\t' + align_len + '\t' + identity + '\t' +
                                             query_start + '\t' + query_end + '\t' + subject_start + '\t' + subject_end
                                             + '\t' + pdb_organism + '\n')
                                g.writelines('$' + query + '\n')
                                g.writelines('%' + subject + '\n')
                                g.writelines('&' + query_ss + '\n')

                                break


def extract_nt_seq(input_dir, tblastn_file, output):
    f = list(open(os.path.join(input_dir, tblastn_file), 'r'))
    g = open(os.path.join(input_dir, output), 'a')

    for i in range(len(f)):
        if f[i][0] == '>':
            acc = [m.split('\t', 2)[1] for m in [f[i]]][0]
            start = int([m.split('\t', 7)[6] for m in [f[i]]][0])
            end = int([m.split('\t', 8)[7] for m in [f[i]]][0])

            Entrez.email = 'michaelwang1335@gmail.com'
            handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text', seq_start=start,
                                   seq_stop=end, api_key='fb157a95528a202add3c042769fb44071807')
            sequence = ''
            handletext = (str(handle.read()))
            splithandle = handletext.split("\n")
            for line in splithandle:
                if line.startswith('>'):
                    pass
                else:
                    sequence += line.strip()

            g.writelines(f[i:i + 4])
            g.writelines('@' + sequence + '\n')

            sleep(0.1)

    g.close()



def correct_ntseq(path, tblastn_ntseq_file, output):
    """parsing tblastn files"""
    f_ntseq = list(open(os.path.join(path, tblastn_ntseq_file), 'r'))
    g_correct_ntseq = open(os.path.join(path, output), 'a')

    for i in range(len(f_ntseq)):
        if f_ntseq[i][0] == '>':
            query = f_ntseq[i + 1][1:-1]
            # Remove '-' in query.
            if '-' in query:
                query = query.replace("-", "")

            subject = f_ntseq[i + 2][1:-1]
            # Remove '-' in subject.
            # Do not remove '*' in subject, will be removed when longest common sequence performed.
            if '-' in subject:
                subject = subject.replace("-", "")

            ss = f_ntseq[i + 3][1:-1]
            ntseq = f_ntseq[i + 4][1:-1]
            # Replace special character in nucleotide sequences to the first choice, based on IUPAC codes.
            IUPAC_dic = {'Y': 'C', 'R': 'A', 'W': 'A', 'S': 'G', 'K': 'T', 'M': 'C', 'D': 'A', 'V': 'A', 'H': 'A',
                         'B': 'C', 'X': 'A', 'N': 'A'}
            ntseq = "".join([IUPAC_dic.get(m, m) for m in ntseq])

            # If 'X's in subject, remove them and corresponding nucleotide sequence.
            if 'X' in subject:
                nt_copy = ntseq
                nt_rev_copy = ntseq

                Xpos = [pos for pos, char in enumerate(subject) if char == 'X']
                Xpos_nt = []
                Xpos_nt_rev = []
                for j in range(len(Xpos)):
                    Xpos_nt.append(Xpos[j] * 3)
                    Xpos_nt.append(Xpos[j] * 3 + 1)
                    Xpos_nt.append(Xpos[j] * 3 + 2)
                    Xpos_nt_rev.append(len(nt_rev_copy) - Xpos[j] * 3 - 1)
                    Xpos_nt_rev.append(len(nt_rev_copy) - Xpos[j] * 3 - 2)
                    Xpos_nt_rev.append(len(nt_rev_copy) - Xpos[j] * 3 - 3)

                subject = subject.replace("X", "")
                nt_copy = ''.join([nt_copy[j] for j in range(len(nt_copy)) if j not in Xpos_nt])
                nt_rev_copy = ''.join([nt_rev_copy[j] for j in range(len(nt_rev_copy)) if j not in Xpos_nt_rev])

                nt_copy_seq = Seq(nt_copy, generic_dna)
                nt_copy_t = nt_copy_seq.translate()

                nt_rev_copy_seq = Seq(nt_rev_copy, generic_dna)
                nt_rev_rev_copy_t = nt_rev_copy_seq.reverse_complement().translate()

                if nt_copy_t == subject:
                    ntseq = str(nt_copy_seq)
                elif nt_rev_rev_copy_t == subject:
                    ntseq = str(nt_rev_copy_seq)
                else:
                    print('Error in: ', f_ntseq[i][1:6], ' after remove "X"')

            ntseq_dna = Seq(ntseq, generic_dna)
            ntseq_dna_rev = ntseq_dna.reverse_complement()

            # If nt sequence read from 3 to 5, mark as '1' for next step use, else mark as '0'.
            if ntseq_dna.translate() == subject:
                g_correct_ntseq.writelines(f_ntseq[i][:-1] + '\t' + '0' + '\n')
                g_correct_ntseq.writelines('$' + query + '\n')
                g_correct_ntseq.writelines('%' + subject + '\n')
                g_correct_ntseq.writelines('&' + ss + '\n')
                g_correct_ntseq.writelines('@' + ntseq_dna + '\n')
            elif ntseq_dna_rev.translate() == subject:
                g_correct_ntseq.writelines(f_ntseq[i][:-1] + '\t' + '1' + '\n')
                g_correct_ntseq.writelines('$' + query + '\n')
                g_correct_ntseq.writelines('%' + subject + '\n')
                g_correct_ntseq.writelines('&' + ss + '\n')
                g_correct_ntseq.writelines('@' + ntseq_dna + '\n')
            # Print error if translated nt sequence not equal to subject sequence.
            else:
                print('Error in: ', f_ntseq[i][1:6])


def longest_substring(s1, s2):
    t = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    l, xl, yl = 0, 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                t[x][y] = t[x - 1][y - 1] + 1
                if t[x][y] > l:
                    l = t[x][y]
                    xl = x
                    yl = y
            else:
                t[x][y] = 0
    # Return sequence length, query_start, query_end, query_common_seq, subject_start, subject_end, subject_common_seq.
    return l, xl - l, xl, s1[xl - l: xl], yl - l, yl, s2[yl - l:yl]



def ntseq_longest_common_string(path, correct_ntseq_file, total_match_ntseq_file, min_len):
    """min_len is the minimum match sequence length."""
    f = list(open(os.path.join(path, correct_ntseq_file), 'r'))
    g = open(os.path.join(path, total_match_ntseq_file), 'a')

    for i in range(len(f)):
        if f[i][0] == '>':
            pdb_id = f[i][1:6]
            acc = [m.split('\t', 2)[1] for m in [f[i]]][0]
            query_start = int([m.split('\t', 5)[4] for m in [f[i]]][0])
            nt_start = int([m.split('\t', 7)[6] for m in [f[i]]][0])
            nt_end = int([m.split('\t', 8)[7] for m in [f[i]]][0])
            organism = [m.split('\t', 9)[8] for m in [f[i]]][0]

            query = f[i + 1][1:-1]
            subject = f[i + 2][1:-1]
            ss = f[i + 3][1:-1]
            ntseq = f[i + 4][1:-1]
            rev_mark = f[i][-2]

            match_len, match_query_start, match_query_end, match_query, match_subject_start, match_subject_end, \
            match_subject = longest_substring(query, subject)

            match_ss = ss[match_query_start:match_query_end]

            if match_len > min_len:
                # if ntseq read from 5 to 3, rev_mark = '0'.
                if rev_mark == '0':
                    match_ntseq_start = match_subject_start * 3
                    match_ntseq_end = match_subject_start * 3 + match_len * 3
                    match_ntseq = ntseq[match_ntseq_start:match_ntseq_end]

                    g.writelines('>' + pdb_id + '\t' + acc + '\t' + str(match_len) + '\t' +
                                 str(query_start + match_query_start) + '\t' +
                                 str(query_start + match_query_start + match_len - 1) + '\t' +
                                 str(nt_start + match_query_start * 3) + '\t' +
                                 str(nt_start + match_query_start * 3 + match_len * 3 - 1) + '\t' + organism + '\t' +
                                 '0' + '\n')
                    g.writelines('$' + match_query + '\n')
                    g.writelines('&' + match_ss + '\n')
                    g.writelines('@' + match_ntseq + '\n')

                # if ntseq read from 3 to 5, rev_mark = '1'.
                elif rev_mark == '1':
                    match_ntseq_end = len(ntseq) - match_subject_start * 3
                    match_ntseq_start = match_ntseq_end - match_len * 3
                    match_ntseq = Seq(ntseq[match_ntseq_start:match_ntseq_end], generic_dna)
                    match_ntSeq = match_ntseq.reverse_complement()

                    # Watch out: careful about the start and end position of nucleotide sequence.
                    g.writelines('>' + pdb_id + '\t' + acc + '\t' + str(match_len) + '\t' +
                                 str(query_start + match_query_start) + '\t' +
                                 str(query_start + match_query_start + match_len - 1) + '\t' +
                                 str(nt_end - match_subject_start * 3 - match_len * 3 + 1) + '\t' +
                                 str(nt_end - match_subject_start * 3) + '\t' + organism + '\t' + '1' + '\n')
                    g.writelines('$' + match_query + '\n')
                    g.writelines('&' + match_ss + '\n')
                    g.writelines('@' + match_ntSeq + '\n')

    g.close()


def aaa_to_a(aaa):
    """Convert 3-digits amino acids to 1-digit."""
    dic = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
           'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
           'TRP': 'W', 'TYR': 'Y'}
    return dic[aaa]


def pdb_rewrite(pdb_match_dir, pdb_match_file, pdb_dir, output_dir):
    """
    Need to download pdb files first.
    Rewrite pdb file, only keep CA and remove alternative points.
    """
    f_pdb_match = list(open(os.path.join(pdb_match_dir, pdb_match_file), 'r'))
    for j in range(len(f_pdb_match)):
        if f_pdb_match[j][0] == '>':
            pdb_id = f_pdb_match[j][1:6]

            f = list(open(os.path.join(pdb_dir, "%s.pdb" % pdb_id[:4]), 'r'))
            g = open(os.path.join(output_dir, "%s_pdb_rewrite" % pdb_id), 'a')

            for i in range(len(f)):
                if f[i][0:4] == 'ATOM' and f[i][21] == pdb_id[4] and f[i][13:15] == 'CA':
                    aaa = f[i][17:20]
                    # SEC, UNK are not consider as amino acid in this research.
                    amino_acid = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN',
                                  'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
                    if aaa in amino_acid:
                        a = aaa_to_a(aaa)
                    else:
                        a = 'X'

                    # 1-digit amino acid + residue# + CA_coordinates.
                    g.writelines(a + '\t' + str(int(f[i][22:26])) + '\t' + str(float(f[i][30:38])) + '\t' +
                                 str(float(f[i][38:46])) + '\t' + str(float(f[i][46:54])) + '\n')
            g.close()


def nt_to_codon(nt):
    """Convert nucleotide codon to special marked amino acid."""
    dic = {'GCT': 'A1', 'GCC': 'A2', 'GCA': 'A3', 'GCG': 'A4', 'TGT': 'C1', 'TGC': 'C2', 'GAT': 'D1', 'GAC': 'D2',
           'GAA': 'E1', 'GAG': 'E2', 'TTT': 'F1', 'TTC': 'F2', 'GGT': 'G1', 'GGC': 'G2', 'GGA': 'G3', 'GGG': 'G4',
           'CAT': 'H1', 'CAC': 'H2', 'ATT': 'I1', 'ATC': 'I2', 'ATA': 'I3', 'AAA': 'K1', 'AAG': 'K2', 'TTA': 'L1',
           'TTG': 'L2', 'CTT': 'L3', 'CTC': 'L4', 'CTA': 'L5', 'CTG': 'L6', 'ATG': 'M1', 'AAT': 'N1', 'AAC': 'N2',
           'CCT': 'P1', 'CCC': 'P2', 'CCA': 'P3', 'CCG': 'P4', 'CAA': 'Q1', 'CAG': 'Q2', 'CGT': 'R1', 'CGC': 'R2',
           'CGA': 'R3', 'CGG': 'R4', 'AGA': 'R5', 'AGG': 'R6', 'TCT': 'S1', 'TCC': 'S2', 'TCA': 'S3', 'TCG': 'S4',
           'AGT': 'S5', 'AGC': 'S6', 'ACT': 'T1', 'ACC': 'T2', 'ACA': 'T3', 'ACG': 'T4', 'GTT': 'V1', 'GTC': 'V2',
           'GTA': 'V3', 'GTG': 'V4', 'TGG': 'W1', 'TAT': 'Y1', 'TAC': 'Y2'}

    # Few ntseq has few special characters, replace to one of the choice.
    special = ['U', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '.', ',']

    if any(alphabet in nt for alphabet in special):
        # Based on IUPAC nucleotide code, replace special character to their first Base choice.
        nt = nt.replace('R', 'A').replace('Y', 'C').replace('S', 'G').replace('S', 'G').replace('W', 'A'). \
            replace('K', 'G').replace('M', 'A').replace('B', 'C').replace('D', 'A').replace('H', 'A'). \
            replace('V', 'A').replace('N', 'A')

    return dic[nt]


def codon_to_pdb_dssp(pdb_rewrite_dir, pdb_rewrite_id, seq_dic, output_dir, bad_dir):
    """
    Convert amino acids into codons in pdb files, add DSSP. seq_dic contains pdb_id and corresponding nt and ss
    sequence.
    """
    f = list(open(os.path.join(pdb_rewrite_dir, "%s_pdb_rewrite" % pdb_rewrite_id), 'r'))
    pdb_seq = ''
    res_num = 0
    res_num_save_line = []
    for i in range(len(f)):
        res_num_next = [j.split('\t', 2)[1] for j in [f[i][:-1]]][0]
        # Watch out: same position mat have more than one residue options: only save the previous one.
        if res_num_next != res_num:
            pdb_seq += f[i][0]
            res_num = res_num_next
            # Save those lines that used for sequence alignment.
            res_num_save_line.append(i)

    ntseq_aa_seq = ''
    ntseq_nt_seq = ''
    ss_seq = ''
    if pdb_rewrite_id.upper() in seq_dic:
        ntseq_aa_seq = seq_dic[pdb_rewrite_id.upper()][0]
        ntseq_nt_seq = seq_dic[pdb_rewrite_id.upper()][2]
        ss_seq = seq_dic[pdb_rewrite_id.upper()][1]

    if ntseq_aa_seq is not '':
        # Find longest common substring.
        match_len, match_query_start, match_query_end, match_query, match_subject_start, match_subject_end, \
        match_subject = longest_substring(pdb_seq, ntseq_aa_seq)

        # Only write to codon_pdb_dssp file if match sequence length > 19.
        if match_len > 19:
            g = open(os.path.join(output_dir, "%s_pdb_codon_dssp" % pdb_rewrite_id), 'a')
            for i in range(match_len):
                nt = ntseq_nt_seq[int(match_subject_start * 3 + i * 3):int(match_subject_start * 3 + i * 3 + 3)]
                ss = ss_seq[int(match_subject_start + i):int(match_subject_start + i + 1)]
                # Convert nucleotide codon to special marked amino acid.
                nt_codon = nt_to_codon(nt)
                g.writelines(nt_codon + '\t' + f[res_num_save_line[match_query_start + i]][:-1] + '\t' + ss + '\n')
            g.close()
        else:
            g_bad = open(os.path.join(bad_dir, 'bad_file'), 'a')
            g_bad.writelines("Match sequence less than 20: " + pdb_rewrite_id + '\n')
    else:
        g_bad = open(os.path.join(bad_dir, 'bad_file'), 'a')
        g_bad.writelines("No ntseq find for: " + pdb_rewrite_id + '\n')


def download_dssp(dssp_url, dssp_id, output_dir):
    """
    If the secondary structure was downloaded from PDB website at the beginning, the alternative way is
    downloading DSSP from ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp. The following function helps to download DSSP
    and add it to corresponding pdb files. DSSP file contains CA coordinates, so no need to download pdb files.
    """
    with closing(request.urlopen("%s/%s.dssp" % (dssp_url, dssp_id))) as r:
        with open(output_dir, 'wb') as f:
            shutil.copyfileobj(r, f)


def codon_to_pdb_dssp_v2(dssp_dir, seq_dic, output_dir, bad_dir):
    """Use dssp files instead of pdb file. Protein sequence, CA and coordinates may slightly different."""
    for key, value in seq_dic.items():
        f = list(open(os.path.join(dssp_dir, "%s_dssp" % key[:4].lower()), 'r'))

        pdb_seq = ''
        res_num_save_line = []
        # Watch out: same position mat have more than one residue options: only save the previous one.
        for i in range(len(f)):
            if f[i][2] == '#' and f[i + 1][11] == key[4]:
                res_num = int(float(f[i + 1][5:10]))
                res_num_save_line.append(i + 1)
                pdb_seq += f[i + 1][13]
                for j in range(i + 2, len(f)):
                    if f[j][11] == key[4]:
                        res_num_next = int(float(f[j][5:10]))
                        if res_num_next != res_num:
                            pdb_seq += f[j][13]
                            res_num = res_num_next
                            # Save those lines that used for sequence alignment.
                            res_num_save_line.append(j)
                break

        ntseq_aa_seq = seq_dic[key][0]
        ntseq_nt_seq = seq_dic[key][2]
        ss_seq = seq_dic[key][1]

        if ntseq_aa_seq is not '':
            # Find longest common substring.
            match_len, match_query_start, match_query_end, match_query, match_subject_start, match_subject_end, \
            match_subject = longest_substring(pdb_seq, ntseq_aa_seq)

            # Only write to codon_pdb_dssp file if match sequence length > 19.
            if match_len > 19:
                g = open(os.path.join(output_dir, "%s_pdb_codon_dssp_v2" % key), 'a')
                for i in range(match_len):
                    nt = ntseq_nt_seq[int(match_subject_start * 3 + i * 3):int(match_subject_start * 3 + i * 3 + 3)]
                    ss = ss_seq[int(match_subject_start + i):int(match_subject_start + i + 1)]
                    # Convert nucleotide codon to special marked amino acid.
                    nt_codon = nt_to_codon(nt)
                    g.writelines(nt_codon + '\t' + f[res_num_save_line[match_query_start + i]][13] + '\t' +
                                 str(int(float(f[res_num_save_line[match_query_start + i]][5:10]))) + '\t' +
                                 str(float(f[res_num_save_line[match_query_start + i]][116:122])) + '\t' +
                                 str(float(f[res_num_save_line[match_query_start + i]][123:129])) + '\t' +
                                 str(float(f[res_num_save_line[match_query_start + i]][130:136])) + '\t' + ss + '\n')
                g.close()
            else:
                g_bad = open(os.path.join(bad_dir, 'bad_file'), 'a')
                g_bad.writelines("Match sequence less than 20: " + key + '\n')
        else:
            g_bad = open(os.path.join(bad_dir, 'bad_file'), 'a')
            g_bad.writelines("No ntseq find for: " + key + '\n')


if __name__ == '__main__':
    """
    First part: Correlate PDB original file to nucleotide sequence.
    """
    ori_dir = 'INDIVIDUAL/PDB_ori'
    tblastn_dir = 'INDIVIDUAL/TBLASTN'
    pdb_ori_id = 'PDBID_pc30_R2.0_Rf0.25_100_1000.txt'
    pdb_ori_csv = 'PDB_pc30_R2.0_Rf0.25_100_1000.csv'

    prepare_pdb_file(ori_dir, pdb_ori_id, pdb_ori_csv, ori_dir, 'PDB_clean.txt')
    tblastn(ori_dir, 'PDB_clean.txt', tblastn_dir)
    extract_correct_tblastn(ori_dir, 'PDB_clean.txt', tblastn_dir, 'PDB_tblastn')
    extract_nt_seq(ori_dir, 'PDB_tblastn', 'PDB_nt')
    correct_ntseq(ori_dir, 'PDB_nt', 'PDB_correct_nt')
    ntseq_longest_common_string(ori_dir, 'PDB_correct_nt', 'PDB_match.txt', 20)

    """
    Second part: Add nucleotide sequence and dssp to pdb files.
    """
    pdb_dir = 'INDIVIDUAL/PDB'
    pdb_rewrite_dir = 'INDIVIDUAL/backup/PDB_rewrite'
    pdb_codon_dssp_dir = 'INDIVIDUAL/PDB_codon_dssp'
    pdb_rewrite(ori_dir, 'PDB_match.txt', pdb_dir, pdb_rewrite_dir)

    ntseq_dic = {}
    f_ntseq = list(open(os.path.join(ori_dir, 'PDB_match.txt'), 'r'))
    for i in range(len(f_ntseq)):
        if f_ntseq[i][0] == '>':
            pdbid = f_ntseq[i][1:6]
            ntseq_dic[pdbid] = []
            ntseq_dic[pdbid].append(f_ntseq[i + 1][1:-1])
            ntseq_dic[pdbid].append(f_ntseq[i + 2][1:-1])
            ntseq_dic[pdbid].append(f_ntseq[i + 3][1:-1])

    pdb_rewrite_list = []
    for filename in glob.glob(os.path.join(pdb_rewrite_dir, "*_pdb_rewrite")):
        pdb_rewrite_list.append(filename[-17:-12])
    for i in range(len(pdb_rewrite_list)):
        codon_to_pdb_dssp(pdb_rewrite_dir, pdb_rewrite_list[i], ntseq_dic, pdb_codon_dssp_dir, ori_dir)
# Python 3.7

"""
File:           Sec_str_analysis
Author:         Shengyuan Wang
Date:           Jun 21, 2019
Last update:    Jun 5, 2020

Purpose:        Codon level protein secondary structure analysis.

Input files:    Protein pdb_codon_dssp files.
Output files:   Secondary structure analysis, including the chi-square test between observed and expected frequency
                for each codon, and create pie plots which show percentage of secondary structure in different
                simplex class.

Notices:
"""

import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import chisquare


def calculate_ss_exp_obs(dssp_path, dssp_files, output_path):
    """
    Calculate the observed and expected frequency of secondary structure for each codon. Secondary structure would be
    classify into 3 (Helix, Sheet, Coil) or 8 (H, I, G, E, B, S, T, C) groups based on DSSP.
    """
    codon_3ss_dic = {}
    codon_8ss_dic = {}

    for dssp_file in dssp_files:
        f_dssp = list(open(os.path.join(dssp_path, dssp_file), 'r'))
        for j in range(len(f_dssp)):
            dssp = f_dssp[j][-2]
            codon = f_dssp[j][:2]
            # Write into 3 dssp classes: helix (G or H or I), strand (E or B), loop (S or T or C).
            if dssp != ' ':
                if codon in codon_3ss_dic:
                    if dssp in ['G', 'H', 'I']:
                        codon_3ss_dic[codon][0] += 1
                    elif dssp in ['E', 'B']:
                        codon_3ss_dic[codon][1] += 1
                    else:
                        codon_3ss_dic[codon][2] += 1

                else:
                    if dssp in ['G', 'H', 'I']:
                        codon_3ss_dic[codon] = [1, 0, 0]
                    elif dssp in ['E', 'B']:
                        codon_3ss_dic[codon] = [0, 1, 0]
                    else:
                        codon_3ss_dic[codon] = [0, 0, 1]
            else:
                if codon in codon_3ss_dic:
                    codon_3ss_dic[codon][2] += 1
                else:
                    codon_3ss_dic[codon] = [0, 0, 1]

            # Write into 8 dssp classes: 310-helix(G), alpha-helix(H), pi-helix(I), beta-bridge(B), beta-bulge(E),
            # turn(T), curvature(S), blank(C).
            dssp_8list = ['G', 'H', 'I', 'B', 'E', 'T', 'S', 'C']
            if dssp != ' ':
                if codon in codon_8ss_dic:
                    if dssp in dssp_8list:
                        codon_8ss_dic[codon][[m for m in range(7) if dssp == dssp_8list[m]][0]] += 1
                    else:
                        codon_8ss_dic[codon][7] += 1
                else:
                    codon_8ss_dic[codon] = [0, 0, 0, 0, 0, 0, 0, 0]
                    codon_8ss_dic[codon][[m for m in range(7) if dssp == dssp_8list[m]][0]] += 1
            else:
                if codon in codon_8ss_dic:
                    codon_8ss_dic[codon][7] += 1
                else:
                    codon_8ss_dic[codon] = [0, 0, 0, 0, 0, 0, 0, 1]

    # Amino acid secondary structure = Sum of all synonymous codons secondary structure.
    aa_3ss_dic = {}
    aa_8ss_dic = {}
    for key in codon_3ss_dic:
        if key[0] not in aa_3ss_dic:
            aa_3ss_dic[key[0]] = codon_3ss_dic[key]
        else:
            aa_3ss_dic[key[0]] = [x + y for x, y in zip(aa_3ss_dic[key[0]], codon_3ss_dic[key])]
    for key in codon_8ss_dic:
        if key[0] not in aa_8ss_dic:
            aa_8ss_dic[key[0]] = codon_8ss_dic[key]
        else:
            aa_8ss_dic[key[0]] = [x + y for x, y in zip(aa_8ss_dic[key[0]], codon_8ss_dic[key])]

    # Expected count = frequency of synonymous codons in amino acid family * sum of each secondary structure in amino
    # acid family.
    codon_3ss_exp_dic = {}
    codon_8ss_exp_dic = {}
    for key in codon_3ss_dic:
        codon_3ss_exp_dic[key] = [sum(codon_3ss_dic[key]) / sum(aa_3ss_dic[key[0]]) * i for i in aa_3ss_dic[key[0]]]
    for key in codon_8ss_dic:
        codon_8ss_exp_dic[key] = [sum(codon_8ss_dic[key]) / sum(aa_8ss_dic[key[0]]) * i for i in aa_8ss_dic[key[0]]]

    # In this version, amino acid expected count = observed count.
    aa_3ss_exp_dic = {}
    aa_8ss_exp_dic = {}
    for key in codon_3ss_exp_dic:
        if key[0] not in aa_3ss_exp_dic:
            aa_3ss_exp_dic[key[0]] = codon_3ss_exp_dic[key]
        else:
            aa_3ss_exp_dic[key[0]] = [x + y for x, y in zip(aa_3ss_exp_dic[key[0]], codon_3ss_exp_dic[key])]
    for key in codon_8ss_exp_dic:
        if key[0] not in aa_8ss_exp_dic:
            aa_8ss_exp_dic[key[0]] = codon_8ss_exp_dic[key]
        else:
            aa_8ss_exp_dic[key[0]] = [x + y for x, y in zip(aa_8ss_exp_dic[key[0]], codon_8ss_exp_dic[key])]

    # Convert dictionary to list for xlsx file writing.
    codon_aa_3ss = []
    obs_alpha = []
    obs_beta = []
    obs_coil = []
    exp_alpha = []
    exp_beta = []
    exp_coil = []
    for key in codon_3ss_dic:
        codon_aa_3ss.append(key)
        obs_alpha.append(codon_3ss_dic[key][0])
        obs_beta.append(codon_3ss_dic[key][1])
        obs_coil.append(codon_3ss_dic[key][2])
        exp_alpha.append(codon_3ss_exp_dic[key][0])
        exp_beta.append(codon_3ss_exp_dic[key][1])
        exp_coil.append(codon_3ss_exp_dic[key][2])
    for key in aa_3ss_dic:
        codon_aa_3ss.append(key)
        obs_alpha.append(aa_3ss_dic[key][0])
        obs_beta.append(aa_3ss_dic[key][1])
        obs_coil.append(aa_3ss_dic[key][2])
        exp_alpha.append(aa_3ss_exp_dic[key][0])
        exp_beta.append(aa_3ss_exp_dic[key][1])
        exp_coil.append(aa_3ss_exp_dic[key][2])

    codon_aa_8ss = []
    obs_8ss = [[], [], [], [], [], [], [], []]
    exp_8ss = [[], [], [], [], [], [], [], []]

    for key in codon_8ss_dic:
        codon_aa_8ss.append(key)
        for i in range(8):
            obs_8ss[i].append(codon_8ss_dic[key][i])
            exp_8ss[i].append(codon_8ss_exp_dic[key][i])
    for key in aa_8ss_dic:
        codon_aa_8ss.append(key)
        for i in range(8):
            obs_8ss[i].append(aa_8ss_dic[key][i])
            exp_8ss[i].append(aa_8ss_exp_dic[key][i])

    df_3ss = pd.DataFrame({'Codon/Amino acid': codon_aa_3ss, 'Obs-Alpha': obs_alpha, 'Obs-Beta': obs_beta,
                           'Obs-Coil': obs_coil, 'Exp-Alpha': exp_alpha, 'Exp-Beta': exp_beta, 'Exp-Coil': exp_coil})
    df_8ss = pd.DataFrame({'Codon/Amino acid': codon_aa_8ss, 'Obs-G': obs_8ss[0], 'Obs-H': obs_8ss[1],
                           'Obs-I': obs_8ss[2], 'Obs-B': obs_8ss[3], 'Obs-E': obs_8ss[4], 'Obs-T': obs_8ss[5],
                           'Obs-S': obs_8ss[6], 'Obs-C': obs_8ss[7], 'Exp-G': exp_8ss[0], 'Exp-H': exp_8ss[1],
                           'Exp-I': exp_8ss[2], 'Exp-B': exp_8ss[3], 'Exp-E': exp_8ss[4], 'Exp-T': exp_8ss[5],
                           'Exp-S': exp_8ss[6], 'Exp-C': exp_8ss[7]})
    with pd.ExcelWriter(os.path.join(output_path, 'codon_aa_2nd_str_exp_obs_id30.xlsx')) as writer:
        df_3ss.to_excel(writer, sheet_name='3ss_obs_exp')
        df_8ss.to_excel(writer, sheet_name='8ss_obs_exp')
    writer.save()


def ss_chi_square(ss_path, ss_file):
    """
    Performe chi-square test between observed and expected frequency.
    """
    xl_3ss = pd.read_excel(os.path.join(ss_path, ss_file), sheet_name='3ss_obs_exp')
    xl_8ss = pd.read_excel(os.path.join(ss_path, ss_file), sheet_name='8ss_obs_exp')

    aa_3ss_chi_square_dic = {}
    for i in range(79):
        if len(xl_3ss['Codon/Amino acid'][i]) == 1:
            obs = [[], [], []]
            exp = [[], [], []]
            aa_3ss_chi_square_dic[xl_3ss['Codon/Amino acid'][i]] = []

            j = i
            while xl_3ss['Codon/Amino acid'][j + 1][0] == xl_3ss['Codon/Amino acid'][i]:
                obs[0].append(xl_3ss['Obs-Alpha'][j + 1])
                obs[1].append(xl_3ss['Obs-Beta'][j + 1])
                obs[2].append(xl_3ss['Obs-Coil'][j + 1])
                exp[0].append(xl_3ss['Exp-Alpha'][j + 1])
                exp[1].append(xl_3ss['Exp-Beta'][j + 1])
                exp[2].append(xl_3ss['Exp-Coil'][j + 1])
                j += 1
                if j > 79:
                    break

            for k in range(3):
                aa_3ss_chi_square_dic[xl_3ss['Codon/Amino acid'][i]].append(float(chisquare(obs[k], f_exp=exp[k])[1]))

    aa_8ss_chi_square_dic = {}
    for i in range(79):
        if len(xl_8ss['Codon/Amino acid'][i]) == 1:
            obs = [[], [], [], [], [], [], [], []]
            exp = [[], [], [], [], [], [], [], []]
            aa_8ss_chi_square_dic[xl_8ss['Codon/Amino acid'][i]] = []

            j = i
            while xl_8ss['Codon/Amino acid'][j + 1][0] == xl_8ss['Codon/Amino acid'][i]:
                sec_list = ['G', 'H', 'I', 'B', 'E', 'T', 'S', 'C']
                for k in range(8):
                    obs[k].append(xl_8ss['Obs-%s' % sec_list[k]][j + 1])
                    exp[k].append(xl_8ss['Exp-%s' % sec_list[k]][j + 1])

                j += 1
                if j > 79:
                    break

            for k in range(8):
                aa_8ss_chi_square_dic[xl_8ss['Codon/Amino acid'][i]].append(float(chisquare(obs[k], f_exp=exp[k])[1]))

    aa3 = []
    pvalue3 = [[], [], []]
    for key in aa_3ss_chi_square_dic:
        aa3.append(key)
        for i in range(3):
            pvalue3[i].append(aa_3ss_chi_square_dic[key][i])

    aa8 = []
    pvalue8 = [[], [], [], [], [], [], [], []]
    for key in aa_8ss_chi_square_dic:
        aa8.append(key)
        for i in range(8):
            pvalue8[i].append(aa_8ss_chi_square_dic[key][i])

    df_3ss = pd.DataFrame({'Amino acid': aa3, 'Alpha': pvalue3[0], 'Beta': pvalue3[1], 'Coil': pvalue3[2]})
    df_8ss = pd.DataFrame({'Amino acid': aa8, 'G': pvalue8[0], 'H': pvalue8[1], 'I': pvalue8[2], 'B': pvalue8[3],
                           'E': pvalue8[4], 'T': pvalue8[5], 'S': pvalue8[6], 'C': pvalue8[7]})
    with pd.ExcelWriter(os.path.join(ss_path, 'ss_aa_chi_square.xlsx')) as writer:
        df_3ss.to_excel(writer, sheet_name='3ss')
        df_8ss.to_excel(writer, sheet_name='8ss')
    writer.save()


def perc_class_2nd(path, file, output_dir):
    """
    Plot the percentage of 3 secondary structure in different simplex classes.
    """
    f = list(open(os.path.join(path, file), 'r'))

    helix = 0
    sheet = 0
    other = 0
    blank = 0
    none = 0
    for i in range(len(f)):
        for j in range(3):
            sec = [x.split('\t', j + 10)[j + 9] for x in [f[i]]][0]
            if sec in ['H', 'G', 'I']:
                helix += 1
            elif sec in ['E', 'B']:
                sheet += 1
            elif sec in ['S', 'T']:
                other += 1
            elif sec == 'M':
                blank += 1
            elif sec == 'N':
                none += 1
        sec = f[i][-2]
        if sec in ['H', 'G', 'I']:
            helix += 1
        elif sec in ['E', 'B']:
            sheet += 1
        elif sec in ['S', 'T']:
            other += 1
        elif sec == 'M':
            blank += 1
        elif sec == 'N':
            none += 1

    labels = ['Helix', 'Sheet', 'Curvature/Turn', 'Blank', 'None']
    sizes = [helix, sheet, other, blank, none]
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'purple']

    plt.pie(sizes, colors=colors, autopct='%1.1f%%', startangle=140)
    plt.axis('equal')
    plt.legend(labels, loc=3)
    fig = plt.gcf()
    fig.savefig(os.path.join(output_dir, 'Pie plot of 3 secondary structure in %s.eps' % file[12:]))
    plt.close(fig)


def perc_class_2nd_8(path, file, output_dir):
    """
    Plot the percentage of 8 secondary structure in different simplex classes.
    """
    f = list(open(os.path.join(path, file), 'r'))

    sec_dic = {'h': 0, 'g': 0, 'i': 0, 'b': 0, 'e': 0, 't': 0, 's': 0, 'c': 0, 'n': 0}
    sec_list = ['H', 'G', 'I', 'B', 'E', 'T', 'S', 'C', 'N']
    for k in range(len(f)):
        for j in range(3):
            sec = [x.split('\t', j + 10)[j + 9] for x in [f[k]]][0]
            for l in sec_list:
                if sec == l:
                    sec_dic[l.lower()] += 1

        sec = f[k][-2]
        for l in sec_list:
            if sec == l:
                sec_dic[l.lower()] += 1

    labels = ['3-10 helix', 'Alpha helix', 'Pi helix', 'Beta bridge', 'Beta bulge', 'Turn', 'Curvature', 'Blank',
              'None']
    sizes = [sec_dic['g'], sec_dic['h'], sec_dic['i'], sec_dic['b'], sec_dic['e'], sec_dic['t'], sec_dic['s'],
             sec_dic['c'], sec_dic['n']]
    label_sorted = [x for _, x in sorted(zip(sizes, labels))]
    sizes_sorted = sorted(sizes)

    plt.pie(sizes_sorted, autopct='%1.1f%%')
    plt.axis('equal')
    plt.legend(label_sorted, loc=3)
    fig = plt.gcf()
    fig.savefig(os.path.join(output_dir, 'Pie plot of 8 secondary structure in %s.eps' % file[12:]))
    plt.close(fig)


if __name__ == '__main__':
    pdb_codon_dssp_dir = 'INDIVIDUAL/PDB_codon_dssp'
    simplex_analysis_dir = 'INDIVIDUAL/simplex_analysis'
    pie_plot_dir = 'INDIVIDUAL/simplex_analysis/Pie plot of sec_str'
    tes_dir = 'INDIVIDUAL/PDB_tessellation'

    # Calculate observed and expected frequency of different secondary structure for each of the codon.
    dssp_files = []
    for filename in glob.glob(os.path.join(pdb_codon_dssp_dir, "*_pdb_codon_dssp")):
        dssp_files.append(filename[71:])
    calculate_ss_exp_obs(pdb_codon_dssp_dir, dssp_files, simplex_analysis_dir)

    # Must sort codon/amino acid based on alphabet in xlsx file before use.
    ss_chi_square(simplex_analysis_dir, 'codon_aa_2nd_str_exp_obs_id30.xlsx')

    # Calculate percentage of secondary structure in different simplex classes.
    sum_simplex = []
    for filename in glob.glob(os.path.join(simplex_analysis_dir, "sum_simplex_class*")):
        perc_class_2nd(simplex_analysis_dir, filename[73:], pie_plot_dir)
        perc_class_2nd_8(simplex_analysis_dir, filename[73:], pie_plot_dir)

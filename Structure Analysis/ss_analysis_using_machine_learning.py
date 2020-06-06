# Python 3.7

"""
File:
Author:         Shengyuan Wang
Date:           Jan 19, 2019
Last update:    Jun 5, 2020

Purpose:        Perform twist beta-strands analysis using machine learning algorithms.

Input files:    Tessellation files,
Output files:   0. Prepare protein summary file
                1. Summarize file of all simplex4 beta-strands and alpha-helix
                2. Summarize file of all 4 consecutive residues in one beta-strands but not in simplex
                3. Terminal/middle decision
                4. Summarize of terminal/middle tessellation statistical analysis
                5. Calculate rotation angle for each 4-residue-beta-strand, plot histogram, perform direction angles
                analysis
                6. Summarize beta-strands variables including edge lengths, average of edge length, tetrahedrality,
                tetrahedron volume, hydropathy, amino acid volume, rotation angle, bend angle, twist angle
                7. Perform chi-square test on synonymous codons to decide the twist beta-strands tendency
                8. Using one-dimentional analysis to distinguish extremely twisted beta-strands
                9. Using K-mean to classify twisted beta-strands and use Silhouette coeeficient to decide the
                separation of different classes.
                10. Perform PCA and SVM to distinguish extremly twisted beta-strands.
                11. Perform Logistic Regression to classify twisted beta-strands.
                12. Perform Decision Tree to clustering beta-sheets based on rotation angle.

Notices:
"""

import os
import glob
import math
import matplotlib.pylab as plt
import pandas as pd

import statistics
import numpy as np
import scipy.stats
from scipy.stats import chisquare
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn import preprocessing
from sklearn.tree import DecisionTreeClassifier  # Import Decision Tree Classifier
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.linear_model import LogisticRegression
from scipy.ndimage.filters import gaussian_filter1d
from sklearn.feature_selection import RFE
from sklearn.utils import resample


def file_summary(pdb_sum_dir, pdb_sum_csv, pdb_codon_dir, pdb_codon_file, beta_dir):
    """
    Extract protein information for each chain.
    """
    pdb_sum = pd.read_csv(os.path.join(pdb_sum_dir, pdb_sum_csv))

    PDB_ID_col = pdb_sum['PDB ID']
    PDB_Chain_col = pdb_sum['Chain ID']
    PDB_Rf_col = pdb_sum['R Free']
    PDB_st_col = pdb_sum['Structure Title']
    PDB_em_col = pdb_sum['Exp. Method']
    PDB_R_col = pdb_sum['Resolution']
    PDB_C_col = pdb_sum['Classification']
    PDB_U_col = pdb_sum['Uniprot Acc']
    PDB_S_col = pdb_sum['Source']
    PDB_EC_col = pdb_sum['EC No']

    sum_PDBID = []
    sum_Rf = []  # R-free
    sum_st = []  # Structure title
    sum_em = []  # Experiment method
    sum_R = []  # Resolution
    sum_C = []  # Classification
    sum_U = []  # Uniprot Accession #
    sum_S = []  # Source
    sum_EC = []  # Enzyme classification #
    sum_CL = []  # Chain length
    sum_seq = []  # Amino acid sequence
    sum_seq_codon = []  # Nucleotide sequence
    sum_ss = []  # Secondary structure

    for i in range(len(PDB_ID_col)):
        pdbid = PDB_ID_col[i] + PDB_Chain_col[i]
        if pdbid in pdb_codon_file:
            sum_PDBID.append(PDB_ID_col[i] + '.' + PDB_Chain_col[i])
            sum_Rf.append(PDB_Rf_col[i])
            sum_st.append(PDB_st_col[i])
            sum_em.append(PDB_em_col[i])
            sum_R.append(PDB_R_col[i])
            sum_C.append(PDB_C_col[i])
            sum_U.append(PDB_U_col[i])
            sum_S.append(PDB_S_col[i])
            sum_EC.append(PDB_EC_col[i])

            f = list(open(os.path.join(pdb_codon_dir, "%s_pdb_codon_dssp" % pdbid), 'r'))
            sum_CL.append(len(f))

            seq = ''
            seq_codon = ''
            ss = ''
            for j in range(len(f)):
                seq += [m.split('\t', 2)[1] for m in [f[j]]][0]
                seq_codon += f[j][:2]
                ss += f[j][-2]

            sum_seq.append(seq)
            sum_seq_codon.append(seq_codon)
            sum_ss.append(ss)

    df = pd.DataFrame({'PDB ID': sum_PDBID, 'Resolution': sum_R, 'R Free': sum_Rf, 'Structure Title': sum_st,
                       'Classification': sum_C, 'EC #': sum_EC, 'Source': sum_S, 'Exp. Method': sum_em,
                       'Uniprot Acc': sum_U, 'Chain Length': sum_CL, 'Amino Acid Sequence': sum_seq,
                       'Nucleotide Sequence': sum_seq_codon, 'Secondary Structure': sum_ss})

    with pd.ExcelWriter(os.path.join(beta_dir, 'protein summary.xlsx')) as writer:
        df.to_excel(writer, sheet_name='Protein Summary', index=False)
    writer.save()


def class4_beta(input_dir, tes_file, output_dir):
    """
    Find all beta-sheets and alpha-helices in simplex class4.
    """
    f = list(open(os.path.join(input_dir, tes_file), 'r'))
    g_alpha = open(os.path.join(output_dir, 'sum_simplex_class4_alpha'), 'a')
    g_beta = open(os.path.join(output_dir, 'sum_simplex_class4_beta'), 'a')

    for i in range(len(f)):
        res_dis = int([m.split('\t', 7)[6] for m in [f[i][:-1]]][0])

        if res_dis == 0:
            dssp_list = []
            for j in range(4):
                dssp_list.append([m.split('\t', j + 16)[j + 15] for m in [f[i][:-1]]][0])

            # 4 residue DSSP must same.
            if len(dssp_list) == 4 and len(set(dssp_list)) == 1:
                # Write PDB_ID, simplex, residue numbers, edge lengths, volume, tetrahedrality into file.
                file_sum = [tes_file[:5], f[i][:8]]

                # Residue numbers
                for j in range(4):
                    file_sum.append([m.split('\t', j + 3)[j + 2] for m in [f[i][:-1]]][0])

                # Edge lengths
                for j in range(6):
                    file_sum.append([m.split('\t', j + 8)[j + 7] for m in [f[i][:-1]]][0])

                file_sum.append([m.split('\t', 14)[13] for m in [f[i][:-1]]][0])  # Volume
                tetrahedrality = [m.split('\t', 15)[14] for m in [f[i][:-1]]][0]

                if 'E' in dssp_list:
                    for k in range(len(file_sum)):
                        g_beta.writelines(file_sum[k] + '\t')
                    g_beta.writelines(tetrahedrality + '\n')

                if 'H' in dssp_list:
                    for k in range(len(file_sum)):
                        g_alpha.writelines(file_sum[k] + '\t')
                    g_alpha.writelines(tetrahedrality + '\n')

    g_alpha.close()
    g_beta.close()


def distance(a, b):
    """Follow functions help to calculate tetrahedron variables with given coordinates."""
    dist = math.sqrt((float(a[0]) - float(b[0])) ** 2 + (float(a[1]) - float(b[1])) ** 2 +
                     (float(a[2]) - float(b[2])) ** 2)
    return str('%.3f' % dist)


def determinant_3x3(m):
    """Calculate the volume of tetrahedron."""
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
            m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1]) +
            m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]))


def subtract(a, b):
    return (float(a[0]) - float(b[0]),
            float(a[1]) - float(b[1]),
            float(a[2]) - float(b[2]))


def tetrahedron_calc_volume(a, b, c, d):
    return str('%.3f' % (abs(determinant_3x3((subtract(a, b), subtract(b, c), subtract(c, d),))) / 6.0))


def tetrahedron_calc_tetrahedrality(l1, l2, l3, l4, l5, l6):
    """Calculate the tetrahedrality of tetrahedron."""
    list = [l1, l2, l3, l4, l5, l6]
    list = [float(i) for i in list]
    average = sum(list) / len(list)
    sum_edge = 0
    for i in range(0, 5):
        for j in range(i + 1, 6):
            sum_edge += (list[i] - list[j]) ** 2
    return str('%.3f' % (sum_edge / ((average ** 2) * 15)))


def tetrahedron_calculation(a, b, c, d):
    """Summary of edges length, volume and tetrahedrality of tetrahedron."""
    l1 = distance(a, b)
    l2 = distance(a, c)
    l3 = distance(a, d)
    l4 = distance(b, c)
    l5 = distance(b, d)
    l6 = distance(c, d)

    volume = tetrahedron_calc_volume(a, b, c, d)
    tetrahedrality = tetrahedron_calc_tetrahedrality(l1, l2, l3, l4, l5, l6)

    return [l1, l2, l3, l4, l5, l6, volume, tetrahedrality]


def find_beta(input_dir, files, output_dir, simplex_file):
    """Find all beta-sheets not observed in simplex."""
    f_simplex = list(open(os.path.join(output_dir, simplex_file), 'r'))
    g = open(os.path.join(output_dir, 'sum_beta_not_in_simplex'), 'a')

    # Save simplex information.
    simplex_dic = {}
    for i in range(len(f_simplex)):
        temp_res_list = []
        for k in range(4):
            temp_res_list.append(int([m.split('\t', k + 3)[k + 2] for m in [f_simplex[i][:-1]]][0]))
        if f_simplex[i][:5] not in simplex_dic:
            simplex_dic[f_simplex[i][:5]] = [sorted(temp_res_list)]
        else:
            simplex_dic[f_simplex[i][:5]].append(sorted(temp_res_list))

    for file in files:
        f = list(open(os.path.join(input_dir, file), 'r'))

        for i in range(len(f) - 3):
            # Check if 4 consecutive residues are all beta-strands.
            if f[i][-2] == f[i + 1][-2] == f[i + 2][-2] == f[i + 3][-2] == 'E':
                res_list = []
                for j in range(4):
                    res_list.append(int([m.split('\t', 3)[2] for m in [f[i + j][:-1]]][0]))

                # If file not in simplex_dic or residues not in simplex_dic[file], write down.
                if (file[:5] not in simplex_dic) or (res_list not in simplex_dic[file[:5]]):
                    cor = [[], [], [], []]
                    for k in range(4):
                        for j in range(3):
                            cor[k].append([m.split('\t', j + 4)[j + 3] for m in [f[i][:-1]]][0])

                    # Calculate edge lengths, volume and tetrahedrality based on 4 consecutive residues.
                    tessellation_list = tetrahedron_calculation(cor[1], cor[2], cor[3], cor[4])

                    sudo_simplex = str(f[i][:2] + f[i + 1][:2] + f[i + 2][:2] + f[i + 3][:2])
                    file_sum = [file[:5], sudo_simplex]
                    for k in range(4):
                        file_sum.append(res_list[k])
                    file_sum += tessellation_list[:-1]

                    for k in range(len(file_sum)):
                        g.writelines(str(file_sum[k]) + '\t')
                    g.writelines(tessellation_list[-1] + '\n')

    g.close()


def find_position(start, end, number):
    """Decide if twisted beta-sheet in the middle or in the terminal of the entire beta-sheet."""
    if number - start <= 2 or end - (number + 3) <= 2:
        return 1
    else:
        return 2


def class4_beta_terminal(input_dir, tes_files, pdb_codon_dssp_dir, output_dir):
    """If terminal/middle need to be decided, use this version."""
    g_terminal = open(os.path.join(output_dir, 'sum_simplex_class4_beta_terminal'), 'a')
    g_middle = open(os.path.join(output_dir, 'sum_simplex_class4_beta_middle'), 'a')

    g_beta = [g_terminal, g_middle]

    for tes_file in tes_files:
        f = list(open(os.path.join(input_dir, tes_file), 'r'))
        for i in range(len(f)):
            res_dis = int([m.split('\t', 7)[6] for m in [f[i][:-1]]][0])

            if res_dis == 0:
                dssp_list = []
                for j in range(4):
                    dssp_list.append([m.split('\t', j + 16)[j + 15] for m in [f[i][:-1]]][0])

                # 4 residue DSSP must same.
                if len(dssp_list) == 4 and len(set(dssp_list)) == 1:
                    # Write PDB_ID, simplex, residue numbers, edge lengths, volume, tetrahedrality into file.
                    file_sum = [tes_file[:5], f[i][:8]]

                    # Residue numbers
                    for j in range(4):
                        file_sum.append([m.split('\t', j + 3)[j + 2] for m in [f[i][:-1]]][0])

                    # Edge lengths
                    for j in range(6):
                        file_sum.append([m.split('\t', j + 8)[j + 7] for m in [f[i][:-1]]][0])

                    file_sum.append([m.split('\t', 14)[13] for m in [f[i][:-1]]][0])  # Volume
                    tetrahedrality = [m.split('\t', 15)[14] for m in [f[i][:-1]]][0]

                if 'E' in dssp_list and len(dssp_list) == 4 and len(set(dssp_list)) == 1:
                    f_pdb = list(open(os.path.join(pdb_codon_dssp_dir, "%s_pdb_codon_dssp" % tes_file[:5]), 'r'))

                    beta_start = -1
                    beta_end = -1
                    for j in range(5, len(f_pdb) - 5):
                        pdb_codon_res = int([m.split('\t', 3)[2] for m in [f_pdb[j][:-1]]][0])
                        if int([m.split('\t', 3)[2] for m in [f[i][:-1]]][0]) == pdb_codon_res:
                            for k in range(1, j + 1):
                                if f_pdb[j - k][-2] == 'E':
                                    pass
                                elif k == j:
                                    beta_start = int([m.split('\t', 3)[2] for m in [f_pdb[j - k][:-1]]][0])
                                    break
                                else:
                                    beta_start = int([m.split('\t', 3)[2] for m in [f_pdb[j - k + 1][:-1]]][0])
                                    break
                            for k in range(1, len(f_pdb) - j):
                                if f_pdb[j + k][-2] == 'E':
                                    pass
                                elif k == len(f_pdb) - j - 1:
                                    beta_end = int([m.split('\t', 3)[2] for m in [f_pdb[j + k][:-1]]][0])
                                    break
                                else:
                                    beta_end = int([m.split('\t', 3)[2] for m in [f_pdb[j + k - 1][:-1]]][0])
                                    break

                    if beta_start != -1 and beta_end != -1:
                        position_num = find_position(beta_start, beta_end,
                                                     int([m.split('\t', 3)[2] for m in [f[i][:-1]]][0]))
                        for k in range(len(file_sum)):
                            g_beta[position_num - 1].writelines(file_sum[k] + '\t')
                        g_beta[position_num - 1].writelines(tetrahedrality + '\n')
    g_terminal.close()
    g_middle.close()


def find_beta_terminal(input_dir, files, output_dir, simplex_file):
    """If terminal/middle need to be decided, use this version."""
    f_simplex = list(open(os.path.join(output_dir, simplex_file), 'r'))
    g_terminal = open(os.path.join(output_dir, 'sum_beta_not_in_simplex_terminal'), 'a')
    g_middle = open(os.path.join(output_dir, 'sum_beta_not_in_simplex_middle'), 'a')
    g = [g_terminal, g_middle]

    # Save simplex information.
    simplex_dic = {}
    for i in range(len(f_simplex)):
        temp_res_list = []
        for k in range(4):
            temp_res_list.append(int([m.split('\t', k + 3)[k + 2] for m in [f_simplex[i][:-1]]][0]))
        if f_simplex[i][:5] not in simplex_dic:
            simplex_dic[f_simplex[i][:5]] = [sorted(temp_res_list)]
        else:
            simplex_dic[f_simplex[i][:5]].append(sorted(temp_res_list))

    for file in files:
        f = list(open(os.path.join(input_dir, file), 'r'))

        for i in range(len(f) - 3):
            # Check if 4 consecutive residues are all beta-strands.
            if f[i][-2] == f[i + 1][-2] == f[i + 2][-2] == f[i + 3][-2] == 'E':
                res_list = []
                for j in range(4):
                    res_list.append(int([m.split('\t', 3)[2] for m in [f[i + j][:-1]]][0]))

                # If file not in simplex_dic or residues not in simplex_dic[file], write down.
                if (file[:5] not in simplex_dic) or (res_list not in simplex_dic[file[:5]]):
                    cor = [[], [], [], []]
                    for k in range(4):
                        for j in range(3):
                            cor[k].append([m.split('\t', j + 4)[j + 3] for m in [f[i][:-1]]][0])

                    # Calculate edge lengths, volume and tetrahedrality based on 4 consecutive residues.
                    tessellation_list = tetrahedron_calculation(cor[1], cor[2], cor[3], cor[4])

                    sudo_simplex = str(f[i][:2] + f[i + 1][:2] + f[i + 2][:2] + f[i + 3][:2])
                    file_sum = [file[:5], sudo_simplex]
                    for k in range(4):
                        file_sum.append(res_list[k])
                    file_sum += tessellation_list[:-1]

                    try:
                        if (f[i - 1][-2] == 'E' and f[i - 2][-2] == 'E' and f[i - 3][-2] == 'E') and \
                                (f[i + 4][-2] == 'E' and f[i + 5][-2] == 'E' and f[i + 6][-2] == 'E'):
                            for k in range(len(file_sum)):
                                g[1].writelines(str(file_sum[k]) + '\t')
                            g[1].writelines(tessellation_list[-1] + '\n')
                        else:
                            for k in range(len(file_sum)):
                                g[0].writelines(str(file_sum[k]) + '\t')
                            g[0].writelines(tessellation_list[-1] + '\n')
                    except IndexError:
                        for k in range(len(file_sum)):
                            g[0].writelines(str(file_sum[k]) + '\t')
                        g[0].writelines(tessellation_list[-1] + '\n')

    g_middle.close()
    g_terminal.close()


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n - 1)
    return m, m - h, m + h


def tessellation_analysis(input_dir, file, filename):
    f = list(open(os.path.join(input_dir, file), 'r'))
    g = open(os.path.join(input_dir, 'sum terminal_middle tessellation analysis'), 'a')

    edge_sum = [0, 0, 0, 0, 0, 0]
    edge_count = 0
    tetra_sum = 0
    volume_sum = 0
    avg_edge = []
    volume_list = []
    tetra_list = []

    for i in range(len(f)):
        edge_list = []
        for j in range(6):
            edge_list.append(float([m.split('\t', j + 7)[j + 6] for m in [f[i][:-1]]][0]))
        sorted_edge_list = sorted(edge_list)
        avg_edge.append(statistics.mean(edge_list))

        for j in range(6):
            edge_sum[j] += sorted_edge_list[j]

        edge_count += 1

        volume = float([m.split('\t', 13)[12] for m in [f[i][:-1]]][0])
        tetrahedrality = float([m.split('\t', 14)[13] for m in [f[i][:-1]]][0])
        volume_list.append(volume)
        tetra_list.append(tetrahedrality)

        volume_sum += volume
        tetra_sum += tetrahedrality

    g.writelines(filename + '\t' + str(edge_count) + '\t')  # Dataset and edge count
    edges = []
    for i in range(6):
        edges.append(edge_sum[i] / edge_count)
        g.writelines('%.3f' % (edge_sum[i] / edge_count) + '\t')  # Average of each edge

    g.writelines('%.3f' % (sum(edge_sum) / (6 * edge_count)) + '\t')                # Total average edge lengths
    g.writelines('%.3f' % statistics.stdev(edges) + '\t')                           # Edge lengths standard deviation
    g.writelines('%.3f' % mean_confidence_interval(avg_edge, 0.95)[1] + '\t')       # Average edge lengths 0.95 CI min
    g.writelines('%.3f' % mean_confidence_interval(avg_edge, 0.95)[2] + '\t')       # Average edge lengths 0.95 CI max
    g.writelines('%.3f' % (volume_sum / edge_count) + '\t')                         # Average volume
    g.writelines('%.3f' % mean_confidence_interval(volume_list, 0.95)[1] + '\t')    # Volume 0.95 CI min
    g.writelines('%.3f' % mean_confidence_interval(volume_list, 0.95)[2] + '\t')    # Volume 0.95 CI max
    g.writelines('%.3f' % (tetra_sum / edge_count) + '\t')                          # Average tetrahedrality
    g.writelines('%.3f' % mean_confidence_interval(tetra_list, 0.95)[1] + '\t')     # Tetrahedrality 0.95 CI min
    g.writelines('%.3f' % mean_confidence_interval(tetra_list, 0.95)[2] + '\n')     # Tetrahedrality 0.95 CI max

    g.close()


def mid_point(c1, c2):
    return [(c1[0] + c2[0]) / 2, (c1[1] + c2[1]) / 2, (c1[2] + c2[2]) / 2]


def plane(c1, c2, c3):
    """Return the coefficient of equation of a plane which determine by three points."""
    p1 = np.array(c1)
    p2 = np.array(c2)
    p3 = np.array(c3)

    # These two vectors are in the plane
    v1 = p3 - p1
    v2 = p2 - p1

    # the cross product is a vector normal to the plane
    cp = np.cross(v1, v2)
    a, b, c = cp

    # This evaluates a * x3 + b * y3 + c * z3 which equals d
    d = -1.0 * np.dot(cp, p3)

    # The equation is ax + by + cz + d = 0.
    return a, b, c, d


def foot(a, b, c, d, point):
    """Return a point of foot on 3D plane."""
    k = (-a * point[0] - b * point[1] - c * point[2] - d) / (a * a + b * b + c * c)
    x2 = a * k + point[0]
    y2 = b * k + point[1]
    z2 = c * k + point[2]

    return [x2, y2, z2]


def findFoot(c1, c2, c3):
    """Return a foot on 2D."""
    co_lambda = ((c1[0] - c2[0]) * (c3[0] - c1[0]) + (c1[1] - c2[1]) * (c3[1] - c1[1]) + (c1[2] - c2[2]) * (
                c3[2] - c1[2])) / ((
                                           c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2)

    x = co_lambda * (c1[0] - c2[0]) + c1[0]
    y = co_lambda * (c1[1] - c2[1]) + c1[1]
    z = co_lambda * (c1[2] - c2[2]) + c1[2]

    return [x, y, z]


def wiki_dihedral(p):
    """
    Formula from Wikipedia article on "Dihedral angle"; formula was removed from the most recent version of article (
    no idea why, the article is a mess at the moment) but the formula can be found in at this permalink to an old
    version of the article:
    https://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=689165217#Angle_between_three_vectors
    uses 1 sqrt, 3 cross products
    """
    p0 = np.array(p[0])
    p1 = np.array(p[1])
    p2 = np.array(p[2])
    p3 = np.array(p[3])

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b0xb1 = np.cross(b0, b1)
    b1xb2 = np.cross(b2, b1)

    b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)

    y = np.dot(b0xb1_x_b1xb2, b1) * (1.0 / np.linalg.norm(b1))
    x = np.dot(b0xb1, b1xb2)

    return np.degrees(np.arctan2(y, x))


def vector_rotation_angle(v1, v2):
    """Angle between two vectors."""
    dot_sign = np.dot(v1, v2)
    cross_sign = np.cross(v1, v2)[0]
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

    if dot_sign > 0 and cross_sign > 0:
        angle = (-1) * np.degrees(angle) + 360
        return angle
    elif dot_sign < 0 and cross_sign > 0:
        angle = 360 - np.degrees(angle)
        return angle
    else:
        return np.degrees(angle)


def vector_bend_angle(v1, v2):
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return np.degrees(angle)


def twist_rotation_bend_angle(c1, c2, c3, c4):
    # c1, c2, c3, c4 are the coordinates list of four consecutive alpha carbon.
    point_dic = dict()
    point_dic['L'] = mid_point(c1, c2)
    point_dic['M'] = mid_point(c2, c3)
    point_dic['N'] = mid_point(c3, c4)
    point_dic['R'] = findFoot(point_dic['M'], point_dic['L'], c2)
    point_dic['P'] = mid_point(point_dic['M'], point_dic['L'])
    point_dic['Q'] = mid_point(point_dic['M'], point_dic['N'])

    vector_dic = dict()
    vector_dic['LM'] = np.array(point_dic['M']) - np.array(point_dic['L'])
    vector_dic['MN'] = np.array(point_dic['N']) - np.array(point_dic['M'])
    vector_dic['Rc2'] = np.array(c2) - np.array(point_dic['R'])

    # The coefficient D of the plane perpendicular to vector LM and point R lies on it.
    plane_dR = vector_dic['LM'][0] * point_dic['R'][0] * (-1) + vector_dic['LM'][1] * point_dic['R'][1] * (-1) + \
               vector_dic['LM'][2] * point_dic['R'][2] * (-1)

    # The plane perpendicular to vector LM and point R lies on it.
    plane_perpendicular_to_LM_R = np.append(vector_dic['LM'], plane_dR)

    # Foot of N on plane which perpendicular to vector LM.
    foot_uR = foot(plane_perpendicular_to_LM_R[0], plane_perpendicular_to_LM_R[1], plane_perpendicular_to_LM_R[2],
                   plane_perpendicular_to_LM_R[3], point_dic['N'])
    point_dic['UR'] = foot_uR

    vector_dic['Ru'] = np.array(np.array(foot_uR) - point_dic['R'])
    vector_dic['LR'] = np.array(point_dic['R']) - np.array(point_dic['L'])
    vector_dic['MR'] = np.array(point_dic['R']) - np.array(point_dic['M'])
    vector_dic['c2M'] = np.array(point_dic['M']) - np.array(c2)

    twist_angle = wiki_dihedral([c2, point_dic['P'], point_dic['Q'], c3])
    bend_angle = vector_bend_angle(vector_dic['LM'], vector_dic['MN'])

    '''
    https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors/16544330#16544330
    Plane embedded in 3D
    One special case is the case where your vectors are not placed arbitrarily, but lie within a plane with a known 
    normal vector n. Then the axis of rotation will be in direction n as well, and the orientation of n will fix an 
    orientation for that axis. In this case, you can adapt the 2D computation above, including n into the determinant 
    to make its size 3×3.

    dot = x1*x2 + y1*y2 + z1*z2
    det = x1*y2*zn + x2*yn*z1 + xn*y1*z2 - z1*y2*xn - z2*yn*x1 - zn*y1*x2
    angle = atan2(det, dot)
    '''

    # Add a reference vector LR to help define the clockwise angle between Rc2 and Ru.
    dot = vector_dic['Rc2'][0] * vector_dic['Ru'][0] + vector_dic['Rc2'][1] * vector_dic['Ru'][1] + vector_dic['Rc2'][
        2] * vector_dic['Ru'][2]
    det = vector_dic['Rc2'][0] * vector_dic['Ru'][1] * vector_dic['LR'][2] + vector_dic['Ru'][0] * vector_dic['LR'][
        1] * vector_dic['Rc2'][2] + vector_dic['LR'][0] * vector_dic['Rc2'][1] * vector_dic['Ru'][2] - \
          vector_dic['Rc2'][
              2] * vector_dic['Ru'][1] * vector_dic['LR'][0] - vector_dic['Ru'][2] * vector_dic['LR'][1] * vector_dic[
              'Rc2'][0] - vector_dic['LR'][2] * vector_dic['Rc2'][1] * vector_dic['Ru'][0]
    angle = math.atan2(det, dot)
    rotation_angle = np.degrees(angle)
    if rotation_angle < 0:
        rotation_angle += 360
    return twist_angle, bend_angle, rotation_angle


def check_in_dic(dic, key):
    if key in dic:
        dic[key] += 1
    else:
        dic[key] = 1


def hist_beta_angles(input_dir, file, pdb_dir, file_name):
    """Plot histogram of each beta angle."""
    f = list(open(os.path.join(input_dir, file), 'r'))
    g_direction = open(os.path.join(input_dir, "direction summary"), 'a')

    twist_list = []
    bend_list = []
    rotation_list = []
    rotate_codon_list = {}
    rotate_edge_list = {}

    for i in range(len(f)):
        id = f[i][:5]
        r1 = [m.split('\t', 3)[2] for m in [f[i]]][0]

        f_pdb = list(open(os.path.join(pdb_dir, "%s_pdb_codon_dssp" % id), 'r'))
        for j in range(len(f_pdb)):
            r_pdb = [m.split('\t', 3)[2] for m in [f_pdb[j]]][0]
            if r_pdb == r1:
                codon = []
                for k in range(4):
                    codon.append(f_pdb[j + k][:2])

                cor1 = []
                cor2 = []
                cor3 = []
                cor4 = []
                for k in range(3):
                    cor1.append(float([m.split('\t', 4 + k)[3 + k] for m in [f_pdb[j]]][0]))
                    cor2.append(float([m.split('\t', 4 + k)[3 + k] for m in [f_pdb[j + 1]]][0]))
                    cor3.append(float([m.split('\t', 4 + k)[3 + k] for m in [f_pdb[j + 2]]][0]))
                    cor4.append(float([m.split('\t', 4 + k)[3 + k] for m in [f_pdb[j + 3]]][0]))

                twist, bend, rotation = twist_rotation_bend_angle(cor1, cor2, cor3, cor4)
                twist_list.append(twist)
                bend_list.append(bend)
                rotation_list.append(rotation)

                for k in range(4):
                    check_in_dic(rotate_codon_list, codon[k])
                for k in range(3):
                    sort_edge = sorted([codon[k], codon[k + 1]])
                    check_in_dic(rotate_edge_list, sort_edge[0] + sort_edge[1])

                break

    codon_frequency_dic = {k: v / total for total in (sum(rotate_codon_list.values()),)
                           for k, v in rotate_codon_list.items()}
    edge_frequency_dic = {k: v / total for total in (sum(rotate_edge_list.values()),)
                          for k, v in rotate_edge_list.items()}

    codon = []
    edge = []
    codon_frequency = []
    edge_frequency = []
    codon_count = []
    edge_count = []
    for key, value in codon_frequency_dic.items():
        codon.append(key)
        codon_frequency.append(value)
        codon_count.append(rotate_codon_list[key])
    for key, value in edge_frequency_dic.items():
        edge.append(key)
        edge_frequency.append(value)
        edge_count.append(rotate_edge_list[key])

    df_codon = pd.DataFrame({'Codon': codon, 'Frequency': codon_frequency, 'Count': codon_count})
    df_edge = pd.DataFrame({'Edge': edge, 'Frequency': edge_frequency, 'Count': edge_count})
    with pd.ExcelWriter(os.path.join(input_dir, 'codon_edge composition %s.xlsx' % file_name)) as writer:
        df_codon.to_excel(writer, sheet_name='Codon', index=False)
        df_edge.to_excel(writer, sheet_name='Edge', index=False)
    writer.save()

    plt.hist2d(rotation_list, bend_list, bins=(100, 100))
    plt.colorbar()
    plt.xlabel('Rotation angle')
    plt.ylabel('Bend angle')
    fig = plt.gcf()
    fig.savefig(os.path.join(input_dir, 'rotation vs bend %s.eps' % file_name))
    plt.close(fig)

    # Twist, bend, rotation angles 3D plot.
    # x = np.array(twist_list)
    # y = np.array(bend_list)
    # z = np.array(rotation_list)
    #
    # xyz = np.vstack([x, y, z])
    # density = stats.gaussian_kde(xyz)(xyz)
    #
    # idx = density.argsort()
    # x, y, z, density = x[idx], y[idx], z[idx], density[idx]*len(twist_list)
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # p = ax.scatter(x, y, z, c=density)
    # fig.colorbar(p)
    # ax.set_xlabel('Twist angle')
    # ax.set_ylabel('Bend angle')
    # ax.set_zlabel('Rotation angle')
    # plt.show()

    right_up = 0
    right_down = 0
    left_down = 0
    left_up = 0
    for i in range(len(rotation_list)):
        direction_angle = rotation_list[i]
        if 90 > direction_angle > 0:
            right_up += 1
        elif 180 > direction_angle > 90:
            right_down += 1
        elif 270 > direction_angle > 180:
            left_down += 1
        elif 360 > direction_angle > 270:
            left_up += 1
    right = right_up + right_down
    left = left_up + left_down
    up = right_up + left_up
    down = left_down + right_down

    g_direction.writelines(file_name + '\t' + str(right_up) + '\t' + str(right_down) + '\t' + str(left_up) + '\t' +
                           str(left_down) + '\t' + str(right) + '\t' + str(left) + '\t' + str(up) + '\t' + str(down)
                           + '\n')
    g_direction.close()


def hist_beta_angles180(input_dir, file, pdb_dir, file_name):
    f = list(open(os.path.join(input_dir, file), 'r'))
    g_direction = open(os.path.join(input_dir, "direction summary"), 'a')
    g_r180 = open(os.path.join(input_dir, "r180 beta"), 'a')

    twist_list = []
    bend_list = []
    rotation_list = []
    rotate_codon_list = {}
    rotate_edge_list = {}

    for i in range(len(f)):
        id = f[i][:5]
        r1 = [m.split('\t', 3)[2] for m in [f[i]]][0]

        f_pdb = list(open(os.path.join(pdb_dir, "%s_pdb_codon_dssp" % id), 'r'))
        for j in range(len(f_pdb)):
            r_pdb = [m.split('\t', 3)[2] for m in [f_pdb[j]]][0]
            if r_pdb == r1:
                codon = []
                for k in range(4):
                    codon.append(f_pdb[j + k][:2])

                cor1 = []
                cor2 = []
                cor3 = []
                cor4 = []
                for k in range(3):
                    cor1.append(float([m.split('\t', k + 4)[k + 3] for m in [f_pdb[j]]][0]))
                    cor2.append(float([m.split('\t', k + 4)[k + 3] for m in [f_pdb[j + 1]]][0]))
                    cor3.append(float([m.split('\t', k + 4)[k + 3] for m in [f_pdb[j + 2]]][0]))
                    cor4.append(float([m.split('\t', k + 4)[k + 3] for m in [f_pdb[j + 3]]][0]))

                twist, bend, rotation = twist_rotation_bend_angle(cor1, cor2, cor3, cor4)

                if rotation <= 180:
                    for k in range(4):
                        check_in_dic(rotate_codon_list, codon[k])
                    for k in range(3):
                        sort_edge = sorted([codon[k], codon[k + 1]])
                        check_in_dic(rotate_edge_list, sort_edge[0] + sort_edge[1])
                    g_r180.writelines(f[i])
                    twist_list.append(twist)
                    bend_list.append(bend)
                    rotation_list.append(rotation)

                break

    codon_frequency_dic = {k: v / total for total in (sum(rotate_codon_list.values()),)
                           for k, v in rotate_codon_list.items()}
    edge_frequency_dic = {k: v / total for total in (sum(rotate_edge_list.values()),)
                          for k, v in rotate_edge_list.items()}

    codon = []
    edge = []
    codon_frequency = []
    edge_frequency = []
    codon_count = []
    edge_count = []
    for key, value in codon_frequency_dic.items():
        codon.append(key)
        codon_frequency.append(value)
        codon_count.append(rotate_codon_list[key])
    for key, value in edge_frequency_dic.items():
        edge.append(key)
        edge_frequency.append(value)
        edge_count.append(rotate_edge_list[key])

    df_codon = pd.DataFrame({'Codon': codon, 'Frequency': codon_frequency, 'Count': codon_count})
    df_edge = pd.DataFrame({'Edge': edge, 'Frequency': edge_frequency, 'Count': edge_count})
    with pd.ExcelWriter(os.path.join(input_dir, 'codon_edge composition %s.xlsx' % file_name)) as writer:
        df_codon.to_excel(writer, sheet_name='Codon', index=False)
        df_edge.to_excel(writer, sheet_name='Edge', index=False)
    writer.save()

    right_up = 0
    right_down = 0
    left_down = 0
    left_up = 0
    for i in range(len(rotation_list)):
        direction_angle = rotation_list[i]
        if 90 > direction_angle > 0:
            right_up += 1
        elif 180 > direction_angle > 90:
            right_down += 1
        elif 270 > direction_angle > 180:
            left_down += 1
        elif 360 > direction_angle > 270:
            left_up += 1
    right = right_up + right_down
    left = left_up + left_down
    up = right_up + left_up
    down = left_down + right_down

    g_direction.writelines(file_name + '\t' + str(right_up) + '\t' + str(right_down) + '\t' + str(left_up) + '\t' +
                           str(left_down) + '\t' + str(right) + '\t' + str(left) + '\t' + str(up) + '\t' + str(down)
                           + '\n')
    g_direction.close()
    g_r180.close()


def simplex_volume_hydropathy(aa):
    aa_volume_dic = {'A': 88.6, 'C': 108.5, 'D': 111.1, 'E': 138.4, 'F': 189.9, 'G': 60.1, 'H': 153.2, 'I': 166.7,
                     'K': 168.6, 'L': 166.7, 'M': 162.9, 'N': 114.1, 'P': 112.7, 'Q': 143.8, 'R': 173.4, 'S': 89,
                     'T': 116.1, 'V': 140, 'W': 227.8, 'Y': 193.6}

    aa_hydropathy_dic = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
                         'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
                         'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3}

    v = []
    h = []
    for i in range(4):
        v.append(aa_volume_dic[aa[i]])
        h.append(aa_hydropathy_dic[aa[i]])
    v_edge = []
    for i in range(3):
        for j in range(i + 1, 4):
            v_edge.append(aa_volume_dic[aa[i]] + aa_volume_dic[aa[j]])

    v_ave = sum(v) / 4
    h_ave = sum(h) / 4
    v_std = statistics.stdev(v)
    h_std = statistics.stdev(h)
    v_max = max(v)
    v_min = min(v)
    aa = [x for _, x in sorted(zip(v, aa))]
    v = sorted(v)
    v_edge = sorted(v_edge)

    return [v_ave, v_std, v_max, v_min, h_ave, h_std], [v[0], v[1], v[2], v[3]], [v_edge[0], v_edge[1], v_edge[2],
                                                                                  v_edge[3], v_edge[4], v_edge[5]], [
               aa[0], aa[1], aa[2], aa[3]]


def beta_angles_files(input_dir, file, pdb_dir, file_name):
    f = list(open(os.path.join(input_dir, file), 'r'))

    twist_list = []
    bend_list = []
    rotation_list = []
    edge_list = [[], [], [], [], [], []]

    ave_list = []
    v_list = []
    t_list = []
    aa_list = [[], [], [], [], [], []]

    v_codon_list = [[], [], [], []]
    v_edge_list = [[], [], [], [], [], []]
    aa_aa_list = [[], [], [], []]

    for i in range(len(f)):
        id = f[i][:5]
        r1 = [m.split('\t', 3)[2] for m in [f[i]]][0]
        simplex = [m.split('\t', 2)[1] for m in [f[i]]][0]
        aa = []
        for j in range(4):
            aa.append(simplex[j * 2])

        f_pdb = list(open(os.path.join(pdb_dir, "%s_pdb_codon_dssp" % id), 'r'))
        for j in range(len(f_pdb)):
            r_pdb = [m.split('\t', 3)[2] for m in [f_pdb[j]]][0]
            if r_pdb == r1:
                cor1 = []
                cor2 = []
                cor3 = []
                cor4 = []
                for k in range(3):
                    cor1.append(float([m.split('\t', 4 + k)[3 + k] for m in [f_pdb[j]]][0]))
                    cor2.append(float([m.split('\t', 4 + k)[3 + k] for m in [f_pdb[j + 1]]][0]))
                    cor3.append(float([m.split('\t', 4 + k)[3 + k] for m in [f_pdb[j + 2]]][0]))
                    cor4.append(float([m.split('\t', 4 + k)[3 + k] for m in [f_pdb[j + 3]]][0]))

                twist, bend, rotation = twist_rotation_bend_angle(cor1, cor2, cor3, cor4)
                twist_list.append(twist)
                bend_list.append(bend)
                rotation_list.append(rotation)

                e = []
                for k in range(6):
                    e.append(float([m.split('\t', k + 7)[k + 6] for m in [f[i]]][0]))
                ave = statistics.mean(e)
                ave_list.append(ave)

                e = sorted(e)

                for k in range(6):
                    edge_list[k].append(e[k])

                v = float([m.split('\t', 13)[12] for m in [f[i]]][0])
                t = float([m.split('\t', 14)[13] for m in [f[i]]][0])
                t_list.append(t)
                v_list.append(v)

                aa_temp, v_temp, v_edge_temp, aa_aa = simplex_volume_hydropathy(aa)
                for k in range(6):
                    aa_list[k].append(aa_temp[k])

                for k in range(4):
                    v_codon_list[k].append(v_temp[k])

                for k in range(6):
                    v_edge_list[k].append(v_edge_temp[k])

                for k in range(4):
                    aa_aa_list[k].append(aa_aa[k])

    df_angles = pd.DataFrame({'E1': edge_list[0], 'E2': edge_list[1], 'E3': edge_list[2], 'E4': edge_list[3],
                              'E5': edge_list[4], 'E6': edge_list[5], 'Ave': ave_list, 'Tetrahedrality': t_list,
                              'Volume': v_list, 'AA_h_ave': aa_list[0], 'AA_h_std': aa_list[1], 'AA_v_ave':
                                  aa_list[2], 'AA_v_std': aa_list[3], 'AA_v_max': aa_list[4], 'AA_v_min':
                                  aa_list[5], 'V0': v_codon_list[0], 'V1': v_codon_list[1], 'V2': v_codon_list[2],
                              'V3': v_codon_list[3], 'V_edge0': v_edge_list[0], 'V_edge1': v_edge_list[1],
                              'V_edge2': v_edge_list[2], 'V_edge3': v_edge_list[3], 'V_edge4': v_edge_list[4],
                              'V_edge5': v_edge_list[5], 'A1': aa_aa_list[0], 'A2': aa_aa_list[1], 'A3': aa_aa_list[2],
                              'A4': aa_aa_list[3], 'Rotation': rotation_list, 'Bend': bend_list, 'Twist': twist_list})
    df_angles.to_csv(os.path.join(input_dir, 'edges_t_v angles %s.csv' % file_name), index=False)


def load_data(path):
    marks_df = pd.read_csv(path, header=0)
    return marks_df


def plot_svc_decision_function(model, ax=None, plot_support=True):
    """Plot the decision function for a 2D SVC"""
    if ax is None:
        ax = plt.gca()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # create grid to evaluate model
    x = np.linspace(xlim[0], xlim[1], 30)
    y = np.linspace(ylim[0], ylim[1], 30)
    Y, X = np.meshgrid(y, x)
    xy = np.vstack([X.ravel(), Y.ravel()]).T
    P = model.decision_function(xy).reshape(X.shape)

    # plot decision boundary and margins
    ax.contour(X, Y, P, colors='k',
               levels=[-1, 0, 1], alpha=0.5,
               linestyles=['--', '-', '--'])

    # plot support vectors
    if plot_support:
        ax.scatter(model.support_vectors_[:, 0],
                   model.support_vectors_[:, 1],
                   s=300, linewidth=1, facecolors='none');
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def chi_square(ss_path, ss_file, sheetname):
    xl_ss_count = pd.read_excel(os.path.join(ss_path, ss_file), sheet_name=sheetname)

    aa_ss_count_chi_square_dic = {}
    for i in [0, 4, 6, 8, 10, 12, 16, 18, 21, 23, 30, 32, 36, 38, 44, 50, 54]:
        aa = xl_ss_count['Codon'][i][0]
        list180 = []
        list0 = []
        for j in range(6):
            if xl_ss_count['Codon'][i + j][0] == aa:
                list180.append(xl_ss_count['180 Count'][i + j])
                list0.append(xl_ss_count['Balance 0 count'][i + j])

        aa_ss_count_chi_square_dic[aa] = float(chisquare(list180, f_exp=list0)[1])
    aa_ss_count_chi_square_dic['Y'] = float(chisquare([xl_ss_count['180 Count'][59], xl_ss_count['180 Count'][60]],
                                                      f_exp=[xl_ss_count['Balance 0 count'][59],
                                                             xl_ss_count['Balance 0 count'][60]])[1])

    aa_list = []
    pvalue = []
    for key, value in aa_ss_count_chi_square_dic.items():
        aa_list.append(key)
        pvalue.append(value)

    df_ss_count = pd.DataFrame({'Amino acid': aa_list, 'p-value': pvalue})
    with pd.ExcelWriter(os.path.join(ss_path, 'syn_codon_twist_beta_pvalue in simplex.xlsx')) as writer:
        df_ss_count.to_excel(writer, sheet_name='pvalue', index=False)
    writer.save()


def bound_acc(data, variable, min_bound, max_bound, coefficient):
    boundary = []
    acc = []
    for j in range(min_bound, max_bound):
        count = 0
        for i in range(len(data[variable])):
            if data[variable][i] > j * coefficient and data['group'][i] == 'Small':
                count += 1
            elif data[variable][i] < j * coefficient and data['group'][i] == 'Big':
                count += 1

        accuracy = count / len(data[variable])
        boundary.append(j * 0.1)
        acc.append(accuracy)

    return boundary, acc


def one_dimension_cluster(input_dir, simplex_file, not_simplex_file, var1_colname, variable1, var2_colname, variable2):
    simplex_data = load_data(os.path.join(input_dir, simplex_file))
    not_simplex_data = load_data(os.path.join(input_dir, not_simplex_file))
    simplex_data['group'] = np.where(simplex_data['Rotation'] > 180, 'Big', simplex_data['Rotation'])
    not_simplex_data['group'] = np.where(not_simplex_data['Rotation'] > 180, 'Big', not_simplex_data['Rotation'])
    simplex_data['group'] = np.where(simplex_data['Rotation'] < 180, 'Small', simplex_data['group'])
    not_simplex_data['group'] = np.where(not_simplex_data['Rotation'] < 180, 'Small', not_simplex_data['group'])

    df1 = simplex_data
    df1_big = df1[df1.group == 'Big']
    df1_small = df1[df1.group == 'Small']

    df1_small_downsampled = resample(df1_small, replace=False, n_samples=len(df1_big), random_state=0)
    df1_downsampled = pd.concat([df1_big, df1_small_downsampled])
    df1_downsampled = df1_downsampled.reset_index()

    df2 = not_simplex_data
    df2_big = df2[df2.group == 'Big']
    df2_small = df2[df2.group == 'Small']

    df2_small_downsampled = resample(df2_small, replace=False, n_samples=len(df2_big), random_state=0)
    df2_downsampled = pd.concat([df2_big, df2_small_downsampled])
    df2_downsampled = df2_downsampled.reset_index()

    # 1D clustering
    boundary, acc = bound_acc(df1_downsampled, var1_colname, 10, 60, 0.1)
    plt.plot(boundary, acc, linewidth=2, label='simplex')

    boundary, acc = bound_acc(df2_downsampled, var1_colname, 10, 60, 0.1)
    plt.plot(boundary, acc, linewidth=2, label='not in simplex')
    plt.legend(loc='best')
    plt.xlabel(variable1)
    plt.ylabel('Accuracy')
    fig = plt.gcf()
    fig.savefig(os.path.join(input_dir, '%s boundary sum2.eps' % variable1))
    plt.close(fig)

    boundary, acc = bound_acc(df1_downsampled, var2_colname, 10, 40, 0.01)
    plt.plot(boundary, acc, linewidth=2, label='simplex')

    boundary, acc = bound_acc(df2_downsampled, var2_colname, 10, 40, 0.01)
    plt.plot(boundary, acc, linewidth=2, label='not in simplex')
    plt.legend(loc='best')
    plt.xlabel(variable2)
    plt.ylabel('Accuracy')
    fig = plt.gcf()
    fig.savefig(os.path.join(input_dir, '%s boundary sum2.eps' % variable2))
    plt.close(fig)


def kmean_silhouette(input_dir, simplex_file, variables, clusters):
    """Clustering beta-sheet rotation angle using Kmean."""
    simplex_data = load_data(os.path.join(input_dir, simplex_file))
    simplex_data['group'] = np.where(simplex_data['Rotation'] > 180, 'Big', simplex_data['Rotation'])
    simplex_data['group'] = np.where(simplex_data['Rotation'] < 180, 'Small', simplex_data['group'])

    # X = feature values, all the columns except the last column
    X = simplex_data.iloc[:, 6:8]  # Change to variable column

    # y = target values, last column of the data frame
    y = simplex_data.iloc[:, -1]

    # Standardizing the features
    x = StandardScaler().fit_transform(X)

    # Perform kmean cluster.
    n_clusters = clusters
    km = KMeans(n_clusters=n_clusters, random_state=1)
    y_pred = km.fit_predict(x)
    y_pred_t = ['Big' if m == 1 else m for m in y_pred]
    y_pred_t = ['Small' if m == 0 else m for m in y_pred_t]
    print(confusion_matrix(y, y_pred_t))
    print(classification_report(y, y_pred_t))

    kmean_ypred = pd.DataFrame({'ypred': y_pred})
    kmeanDf = pd.concat([X, kmean_ypred], axis=1)

    plt.scatter(
        kmeanDf.loc[kmeanDf['ypred'] == 1, 'Ave'], kmeanDf.loc[kmeanDf['ypred'] == 1, 'Tetrahedrality'],
        s=1, c='red',
        marker='o',
        label='cluster 1'
    )

    plt.scatter(
        kmeanDf.loc[kmeanDf['ypred'] == 0, 'Ave'], kmeanDf.loc[kmeanDf['ypred'] == 0, 'Tetrahedrality'],
        s=1, c='green',
        marker='o',
        label='cluster 2'
    )
    plt.xlabel('Average of edge lengths')
    plt.ylabel('Tetrahedrality')
    fig = plt.gcf()
    fig.savefig(os.path.join(input_dir, '%s kmean predict.eps' % variables))
    plt.close(fig)

    plt.scatter(
        kmeanDf.loc[simplex_data['Rotation'] > 180, 'Ave'], kmeanDf.loc[simplex_data['Rotation'] > 180, 'Tetrahedrality'],
        s=1, c='red',
        marker='o',
        label='cluster 1'
    )

    plt.scatter(
        kmeanDf.loc[simplex_data['Rotation'] < 180, 'Ave'], kmeanDf.loc[simplex_data['Rotation'] < 180, 'Tetrahedrality'],
        s=1, c='green',
        marker='o',
        label='cluster 2'
    )
    plt.xlabel('Average of edge lengths')
    plt.ylabel('Tetrahedrality')
    fig = plt.gcf()
    fig.savefig(os.path.join(input_dir, '%s kmean actual.eps' % variables))
    plt.close(fig)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed clusters
    silhouette_avg = silhouette_score(x, y_pred)

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(x, y_pred)

    fig, ax = plt.subplots()

    ax.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax.set_ylim([0, len(x) + (n_clusters + 1) * 10])

    y_lower = 10
    for i, j, k in [0, 'Big', 'green'], [1, 'Small', 'red']:
        # Aggregate the silhouette scores for samples belonging to cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[y_pred == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = k
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.095, y_lower + 0.5 * size_cluster_i, j)

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title("The average of silhouette_score is %s." % str('%.3f' % silhouette_avg))
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="black", linestyle="--")

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    fig.savefig(os.path.join(input_dir, 'Silhouette plot.eps'))
    plt.close(fig)


def PCA_SVM(input_dir, simplex_file, variables, svm_kernal):
    """
    Exploratory data analysis using PCA.
    Classify beta-sheets using SVM.
    """
    simplex_data = load_data(os.path.join(input_dir, simplex_file))
    simplex_data['group'] = np.where(simplex_data['Rotation'] > 180, 'Big', simplex_data['Rotation'])
    simplex_data['group'] = np.where(simplex_data['Rotation'] < 180, 'Small', simplex_data['group'])

    # X = feature values, all the columns except the last column
    X = simplex_data.iloc[:, 6:8]  # Change to variable column

    # Standardizing the features
    x = StandardScaler().fit_transform(X)

    # Perform PCA analysis.
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])
    finalDf = pd.concat([principalDf, simplex_data['group']], axis=1)

    # Perfrom SVM on PCA results.
    X_train, X_test, y_train, y_test = train_test_split(principalDf, simplex_data['group'], test_size=0.20,
                                                        random_state=0)
    svclassifier = SVC(kernel=svm_kernal)
    svclassifier.fit(X_train, y_train)
    y_pred = svclassifier.predict(X_test)
    print(confusion_matrix(y_test, y_pred))
    print(classification_report(y_test, y_pred))

    markers = ['o', 'x']
    targets = ['Big', 'Small']
    for target, m in zip(targets, markers):
        indicesToKeep = finalDf['group'] == target
        plt.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                    , finalDf.loc[indicesToKeep, 'principal component 2']
                    , marker=m
                    , s=10)
    plot_svc_decision_function(svclassifier)
    plt.scatter(svclassifier.support_vectors_[:, 0], svclassifier.support_vectors_[:, 1],
                s=300, lw=1, facecolors='none')
    plt.xlabel('Principal component 1')
    plt.ylabel('Principal component 2')
    plt.legend(['Rotation angle > 180', 'Rotation angle < 180'])
    fig = plt.gcf()
    fig.savefig(os.path.join(input_dir, '%s PCA_SVM %s.eps' % (variables, svm_kernal)))
    plt.close(fig)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.set_title('2 component PCA', fontsize=20)
    targets = ['Big', 'Small']
    markers = ['o', 'x']
    for target, m in zip(targets, markers):
        indicesToKeep = finalDf['group'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                   , finalDf.loc[indicesToKeep, 'principal component 2']
                   , marker=m
                   , s=10)
    ax.legend(['Rotation angle > 180', 'Rotation angle < 180'])
    ax.grid()
    fig = plt.gcf()
    fig.savefig(os.path.join(input_dir, '%s PCA.eps' % variables))
    plt.close(fig)


def logistic_regression(input_dir, csv_file, output_name):
    """
    Classify beta-sheets using Logistic Regression.
    """
    data = pd.read_csv(os.path.join(input_dir, csv_file))

    data['group'] = np.where(data['Rotation'] > 180, 0, data['Rotation'])
    data['group'] = np.where(data['Rotation'] < 180, 1, data['group'])
    y = data.group

    # columns = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'Ave', 'Tetrahedrality', 'Volume', 'AA_h_ave', 'AA_h_std',
    #               'AA_v_ave', 'AA_v_std', 'AA_v_max', 'AA_v_min', 'V0', 'V1', 'V2', 'V3', 'V_edge0', 'V_edge1',
    #               'V_edge2', 'V_edge3', 'V_edge4', 'V_edge5']

    data1 = data[['E3', 'E6']]
    data2 = data[['A1', 'A2', 'A3', 'A4']]

    le = preprocessing.LabelEncoder()
    data2 = data2.apply(le.fit_transform)

    x = np.concatenate((data1, data2), axis=1)

    X_train, X_test, y_train, y_test = train_test_split(data1, y, test_size=0.2, random_state=0)

    logreg = LogisticRegression(solver='sag', max_iter=1000)

    # Feature selection
    rfe = RFE(logreg, 2)
    rfe = rfe.fit(X_train, y_train.values.ravel())
    print(rfe.support_)
    print(rfe.ranking_)

    logreg.fit(X_train, y_train.values.ravel())
    y_pred = logreg.predict(X_test)

    print('Accuracy of logistic regression classifier on test set: {:.2f}'.format(logreg.score(X_test, y_test)))
    print(confusion_matrix(y_test, y_pred))
    print(classification_report(y_test, y_pred))

    logit_roc_auc = roc_auc_score(y_test, logreg.predict(X_test))
    fpr, tpr, thresholds = roc_curve(y_test, logreg.predict_proba(X_test)[:, 1])

    tpr_smoothed = gaussian_filter1d(tpr, sigma=2)

    X2_train, X2_test, y2_train, y2_test = train_test_split(X2, y, test_size=0.2, random_state=0)
    logreg = LogisticRegression(solver='sag', max_iter=1000)
    logreg.fit(X2_train, y2_train.values.ravel())
    y2_pred = logreg.predict(X2_test)

    logit_roc_auc2 = roc_auc_score(y2_test, logreg.predict(X2_test))
    fpr2, tpr2, thresholds2 = roc_curve(y2_test, logreg.predict_proba(X2_test)[:, 1])
    tpr2_smoothed = gaussian_filter1d(tpr2, sigma=2)

    plt.figure()
    plt.plot(fpr2, tpr2_smoothed, label='All 25 tessellation parameters (area = %0.2f)' % logit_roc_auc2)
    plt.plot(fpr, tpr_smoothed, label='Average of edge lengths and tetrahedrality (area = %0.2f)' % logit_roc_auc)
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(input_dir, output_name))


def decision_tree(input_dir, csv_file):
    """
    Clustering beta-sheets using Decision Tree based on rotation angle.
    """
    data = pd.read_csv(os.path.join(input_dir, csv_file))

    data['group'] = np.where(data['Rotation'] > 180, 0, data['Rotation'])
    data['group'] = np.where(data['Rotation'] < 180, 1, data['group'])

    data_big = data[data.group == 0]
    data_small = data[data.group == 1]

    data_small_downsampled = resample(data_small, replace=False, n_samples=len(data_big), random_state=0)
    data_downsampled = pd.concat([data_big, data_small_downsampled])
    data = data_downsampled.reset_index()
    y = data.group

    columns = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'Ave', 'Tetrahedrality', 'Volume', 'AA_h_ave', 'AA_h_std',
               'AA_v_ave', 'AA_v_std', 'AA_v_max', 'AA_v_min', 'V0', 'V1', 'V2', 'V3', 'V_edge0', 'V_edge1',
               'V_edge2', 'V_edge3', 'V_edge4', 'V_edge5']
    data1 = data[['Ave', 'Tetrahedrality']]
    data2 = data[['A1', 'A2', 'A3', 'A4']]

    le = preprocessing.LabelEncoder()
    data2 = data2.apply(le.fit_transform)

    x = np.concatenate((data1, data2), axis=1)

    sum_accuracy = 0
    for i in range(5):
        # Split dataset into training set and test set
        X_train, X_test, y_train, y_test = train_test_split(data2, y, test_size=0.2, random_state=i)

        # Create Decision Tree classifer object
        clf = DecisionTreeClassifier()

        # Train Decision Tree Classifer
        clf = clf.fit(X_train, y_train)

        # Predict the response for test dataset
        y_pred = clf.predict(X_test)

        accuracy = metrics.accuracy_score(y_test, y_pred)
        sum_accuracy += accuracy

    avg_accuracy = sum_accuracy / 5

    print(avg_accuracy)


if __name__ == '__main__':
    pdb_ori_dir = 'INDIVIDUAL/PDB_ori'
    tes_dir = 'INDIVIDUAL/PDB_tessellation'
    pdb_codon_dssp_dir = 'INDIVIDUAL/PDB_codon_dssp'
    beta_dir = 'INDIVIDUAL/TWIST_beta_strand'

    # Prepare protein summary file.
    pdb_codon_files = []
    for filename in glob.glob(os.path.join(pdb_codon_dssp_dir, "*_pdb_codon_dssp")):
        pdb_codon_files.append(filename[-20:-15])

    file_summary(pdb_ori_dir, 'PDB_pc30_R2.0_Rf0.25_100_1000.csv', pdb_codon_dssp_dir, pdb_codon_files, beta_dir)

    # Find all beta-strand and alpha-helix which formed by 4 consecutive residues and tessellated in one simplex.
    for filename in glob.glob(os.path.join(tes_dir, "*_tessellation")):
        class4_beta(tes_dir, filename[-18:], beta_dir)

    # Find those beta-strands formed by 4 consecutive residues but not in simplex.
    pdb_codon_dssp_files = []
    for filename in glob.glob(os.path.join(pdb_codon_dssp_dir, "*_pdb_codon_dssp")):
        pdb_codon_dssp_files.append(filename[-20:])
    find_beta(pdb_codon_dssp_dir, pdb_codon_dssp_files, beta_dir, 'sum_simplex_class4_alpha')

    # Classify beta-strand into middle or terminal.
    tes_files = []
    for filename in glob.glob(os.path.join(tes_dir, "*_tessellation")):
        tes_files.append(filename[-18:])
    class4_beta_terminal(tes_dir, tes_files, pdb_codon_dssp_dir, beta_dir)

    # Classify those beta-strands formed by 4 consecutive residues but not in simplex into terminal or middle.
    pdb_codon_dssp_files = []
    for filename in glob.glob(os.path.join(pdb_codon_dssp_dir, "*_pdb_codon_dssp")):
        pdb_codon_dssp_files.append(filename[-20:])
    find_beta_terminal(pdb_codon_dssp_dir, pdb_codon_dssp_files, beta_dir, 'sum_simplex_class4_beta')

    # Statistical analysis on terminal/middle datasets.
    tessellation_analysis(beta_dir, 'sum_simplex_class4_beta_terminal', 'Sp_terminal')
    tessellation_analysis(beta_dir, 'sum_simplex_class4_beta_middle', 'Sp_middle')
    tessellation_analysis(beta_dir, 'sum_beta_not_in_simplex_terminal', 'Nm_terminal')
    tessellation_analysis(beta_dir, 'sum_beta_not_in_simplex_middle', 'Nm_middle')
    tessellation_analysis(beta_dir, 'sum_simplex_class4_alpha', 'Sp_alpha')
    tessellation_analysis(beta_dir, 'sum_alpha_not_in_simplex', 'Nm_alpha')

    # Classify based on rotation angle less or greater than 180.
    hist_beta_angles(beta_dir, 'sum_beta_not_in_simplex', pdb_codon_dssp_dir, 'not in simplex')
    hist_beta_angles(beta_dir, 'sum_simplex_class4_beta', pdb_codon_dssp_dir, 'within simplex')
    hist_beta_angles180(beta_dir, 'sum_beta_all', pdb_codon_dssp_dir, 'greater r180')
    hist_beta_angles180(beta_dir, 'sum_simplex_class4_beta', pdb_codon_dssp_dir, 'greater r180 in simplex')
    hist_beta_angles180(beta_dir, 'sum_simplex_class4_beta', pdb_codon_dssp_dir, 'less r180 in simplex')

    # Create csv file with each beta-strand information.
    beta_angles_files(beta_dir, 'sum_simplex_class4_beta', pdb_codon_dssp_dir, 'within simplex')
    beta_angles_files(beta_dir, 'sum_beta_not_in_simplex', pdb_codon_dssp_dir, 'not in simplex')

    # Use average of edge lengths and tetrahedrality to classify data.
    one_dimension_cluster(beta_dir, 'edges_t_v angles within simplex.csv', 'edges_t_v angles not in simplex.csv',
                          'Ave', 'Average of edge lengths', 'Tetrahedrality', 'Tetrahedrality')

    # Perform chi-square test between greater than 180 count and less than 180 count on synonymous codons.
    chi_square(beta_dir, 'codon_edge composition greater r180 in simplex.xlsx', 'Comparison')

    # Perform Kmean on selected variables.
    kmean_silhouette(beta_dir, 'edges_t_v angles within simplex.csv', 'Ave_Tetrahedra', 2)

    # Perform PCA on selected variables, use SVM to find two clusters.
    PCA_SVM(beta_dir, 'edges_t_v angles within simplex.csv', 'Ave_Tetrahedra', 'rbf')

    logistic_regression(beta_dir, 'edges_t_v angles within simplex.csv', 'LR_ROC')
    decision_tree(beta_dir, 'edges_t_v angles within simplex.csv')

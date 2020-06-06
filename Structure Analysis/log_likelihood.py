# Python 3.7

"""
File:           log_likelihood.py
Author:         Shengyuan Wang
Date:           Jun 16, 2019
Last update:    Jun 4, 2020

Purpose:        calculate the log-likelihood of each simplex.

Input files:    Tessellation files.
Output files:   Log-likelihood of each kind of simplex.

Notices:        1. log-likelihood of ijkl = f(ijkl) / p(ijkl), where p(ijkl) = ca(i)a(j)a(k)a(l), where a(m) is the
                frequency of m in the training set.
                2. This version1 calculate the log-likelihood of simplices by using summarized protein files (i.e. all
                proteins in 1 training set) and calculate the frequencies of codons in simplices (i.e. do not remove
                the repeats).
"""


import os
import glob
import collections
import math
import statistics
import pandas as pd


def write_summary_simplex(tes_dir, tes_files, llh_dir, cutoff):
    """Combine sum_simplex_class files, add class in the end of each line."""
    g = open(os.path.join(llh_dir, 'sum_all_tessellation_cutoff_%s' % str(cutoff)), 'a')

    for i in range(len(tes_files)):
        f = list(open(os.path.join(tes_dir, tes_files[i]), 'r'))

        for j in range(len(f)):
            edge_list = []
            for k in range(6):
                edge = float([m.split('\t', k+8)[k+7] for m in [f[j]]][0])
                edge_list.append(edge)

            if all(x < cutoff for x in edge_list):
                g.writelines(tes_files[i][:5] + '\t' + f[j])
    g.close()


def check_key_in_dic(dic, key_list):
    for i in range(len(key_list)):
        if key_list[i] in dic:
            dic[key_list[i]] += 1
        else:
            dic[key_list[i]] = 1


def simplex_permutation_factor(one_simplex):
    """Calculate permutation factor of each simplex."""
    one_simplex1 = [one_simplex[0:2], one_simplex[2:4], one_simplex[4:6], one_simplex[6:8]]
    counter = collections.Counter(one_simplex1)
    if len(counter) == 1:
        return 1
    elif len(counter) == 2:
        if counter[one_simplex[0:2]] == 2:
            return 6
        else:
            return 4
    elif len(counter) == 3:
        return 12
    elif len(counter) == 4:
        return 24


def simplex_permutation_factor_aa(one_simplex):
    one_simplex1 = [one_simplex[0], one_simplex[1], one_simplex[2], one_simplex[3]]
    counter = collections.Counter(one_simplex1)
    if len(counter) == 1:
        return 1
    elif len(counter) == 2:
        if counter[one_simplex[0:2]] == 2:
            return 6
        else:
            return 4
    elif len(counter) == 3:
        return 12
    elif len(counter) == 4:
        return 24


def calculate_simplex_log_likelihood(llh_dir, sum_simplex_file):
    """Calculate simplex log-likelihood on codon level."""
    f_sum_simplex = list(open(os.path.join(llh_dir, sum_simplex_file), 'r'))

    codon_dic = {}
    simplex_dic = {}
    simplex_v_dic = {}
    simplex_t_dic = {}

    for i in range(len(f_sum_simplex)):
        simplex = [m.split('\t', 3)[2] for m in [f_sum_simplex[i]]][0]
        codon1 = simplex[:2]
        codon2 = simplex[2:4]
        codon3 = simplex[4:6]
        codon4 = simplex[6:8]
        volume = float([m.split('\t', 15)[14] for m in [f_sum_simplex[i]]][0])
        tetrahedrality = float([m.split('\t', 16)[15] for m in [f_sum_simplex[i]]][0])

        # Add 1 to value if codon or simplex appears.
        check_key_in_dic(codon_dic, [codon1, codon2, codon3, codon4])
        check_key_in_dic(simplex_dic, [simplex])

        # Create volume and tetrahedrality dictionary for simplices.
        if simplex not in simplex_v_dic:
            simplex_v_dic[simplex] = [volume]
        else:
            simplex_v_dic[simplex].append(volume)

        if simplex not in simplex_t_dic:
            simplex_t_dic[simplex] = [tetrahedrality]
        else:
            simplex_t_dic[simplex].append(tetrahedrality)

    # Calculate the frequency of each codon, simplex.
    codon_frequency_dic = {k: v / total for total in (sum(codon_dic.values()),) for k, v in codon_dic.items()}
    simplex_frequency_dic = {k: v / total for total in (sum(simplex_dic.values()),) for k, v in simplex_dic.items()}

    # Calculate log-likelihood for each simplex.
    simplex_llh_dic = {}
    for key in simplex_frequency_dic:
        simplex_llh_dic[key] = math.log10(simplex_frequency_dic[key] / (simplex_permutation_factor(key) *
                                codon_frequency_dic[key[:2]] * codon_frequency_dic[key[2:4]] *
                                codon_frequency_dic[key[4:6]] * codon_frequency_dic[key[6:8]]))

    # Calculate the average volume and tetrahedrality of each simplex.
    simplex_average_v_dic = {}
    for key, value in simplex_v_dic.items():
        simplex_average_v_dic[key] = statistics.mean(value)

    simplex_average_t_dic = {}
    for key, value in simplex_t_dic.items():
        simplex_average_t_dic[key] = statistics.mean(value)

    # Write into xlsx file with Simplex and Codon.
    simplex = []
    simplex_llh = []
    simplex_frequency = []
    simplex_count = []
    simplex_average_v = []
    simplex_average_t = []
    codon = []
    codon_count = []
    codon_frequency = []
    for key in simplex_llh_dic:
        simplex.append(key)
        simplex_llh.append(simplex_llh_dic[key])
        simplex_frequency.append(simplex_frequency_dic[key])
        simplex_count.append(simplex_dic[key])
        simplex_average_v.append(simplex_average_v_dic[key])
        simplex_average_t.append(simplex_average_t_dic[key])
    for key in codon_frequency_dic:
        codon.append(key)
        codon_count.append(codon_dic[key])
        codon_frequency.append(codon_frequency_dic[key])

    # Excel cannot handle files with more than 1M rows, either split into two xlsx files, or wirte into txt file.
    # df_simplex = pd.DataFrame({'Simplex': simplex, 'Log-likelihood': simplex_llh, 'Frequency': simplex_frequency,
    #                            'Count': simplex_count, 'Average volume': simplex_average_v, 'Average tetrahedrality':
    #                                simplex_average_t})
    # df_codon = pd.DataFrame({'Codon': codon, 'Frequency': codon_frequency, 'Count': codon_count})
    # with pd.ExcelWriter(os.path.join(llh_dir, 'simplex_llh_cutoff_12.xlsx')) as writer:
    #     df_simplex.to_excel(writer, sheet_name='Simplex')
    #     df_codon.to_excel(writer, sheet_name='Codon')
    # writer.save()

    # g = open(os.path.join(llh_dir, 'simplex_llh_cutoff8.txt'), 'a')
    # for i in range(len(simplex)):
    #     g.writelines(simplex[i] + '\t' + str("%.7f" % simplex_llh[i]) + '\t' + str("%.7f" % simplex_frequency[i]) +
    #                  '\t' + str(simplex_count[i]) + '\t' + str("%.3f" % simplex_average_v[i]) + '\t' +
    #                  str("%.3f" % simplex_average_t[i]) + '\n')
    # g.close()

    g = open(os.path.join(llh_dir, 'codon_frequency_cutoff8.txt'), 'a')
    for i in range(len(codon)):
        g.writelines(codon[i] + '\t' + str("%.4f" % codon_frequency[i]) + '\t' + str(codon_count[i]) + '\n')
    g.close()


def calculate_simplex_log_likelihood_aa(llh_dir, sum_simplex_file):
    """Calculate simplex log-likelihood on amino acid level."""
    f_sum_simplex = list(open(os.path.join(llh_dir, sum_simplex_file), 'r'))

    aa_dic = {}
    simplex_dic = {}
    simplex_v_dic = {}
    simplex_t_dic = {}

    for i in range(len(f_sum_simplex)):
        simplex = [m.split('\t', 3)[2] for m in [f_sum_simplex[i]]][0]
        simplex_aa = str(simplex[0]+simplex[2]+simplex[4]+simplex[6])
        aa1 = simplex[0]
        aa2 = simplex[2]
        aa3 = simplex[4]
        aa4 = simplex[6]
        volume = float([m.split('\t', 15)[14] for m in [f_sum_simplex[i]]][0])
        tetrahedrality = float([m.split('\t', 16)[15] for m in [f_sum_simplex[i]]][0])

        # Add 1 to value if codon or simplex appears.
        check_key_in_dic(aa_dic, [aa1, aa2, aa3, aa4])
        check_key_in_dic(simplex_dic, [simplex_aa])

        # Create volume and tetrahedrality dictionary for simplices.
        if simplex_aa not in simplex_v_dic:
            simplex_v_dic[simplex_aa] = [volume]
        else:
            simplex_v_dic[simplex_aa].append(volume)

        if simplex_aa not in simplex_t_dic:
            simplex_t_dic[simplex_aa] = [tetrahedrality]
        else:
            simplex_t_dic[simplex_aa].append(tetrahedrality)

    # Calculate the frequency of each aa, simplex.
    aa_frequency_dic = {k: v / total for total in (sum(aa_dic.values()),) for k, v in aa_dic.items()}
    simplex_frequency_dic = {k: v / total for total in (sum(simplex_dic.values()),) for k, v in simplex_dic.items()}

    # Calculate log-likelihood for each simplex.
    simplex_llh_dic = {}
    for key in simplex_frequency_dic:
        simplex_llh_dic[key] = math.log10(simplex_frequency_dic[key] / (simplex_permutation_factor_aa(key) *
                                aa_frequency_dic[key[0]] * aa_frequency_dic[key[1]] *
                                aa_frequency_dic[key[2]] * aa_frequency_dic[key[3]]))

    # Calculate the average volume and tetrahedrality of each simplex.
    simplex_average_v_dic = {}
    for key, value in simplex_v_dic.items():
        simplex_average_v_dic[key] = statistics.mean(value)

    simplex_average_t_dic = {}
    for key, value in simplex_t_dic.items():
        simplex_average_t_dic[key] = statistics.mean(value)

    # Write into xlsx file with Simplex and aa.
    simplex = []
    simplex_llh = []
    simplex_frequency = []
    simplex_count = []
    simplex_average_v = []
    simplex_average_t = []
    aa = []
    aa_count = []
    aa_frequency = []
    for key in simplex_llh_dic:
        simplex.append(key)
        simplex_llh.append(simplex_llh_dic[key])
        simplex_frequency.append(simplex_frequency_dic[key])
        simplex_count.append(simplex_dic[key])
        simplex_average_v.append(simplex_average_v_dic[key])
        simplex_average_t.append(simplex_average_t_dic[key])
    for key in aa_frequency_dic:
        aa.append(key)
        aa_count.append(aa_dic[key])
        aa_frequency.append(aa_frequency_dic[key])

    # Excel cannot handle files with more than 1M rows, either split into two xlsx files, or wirte into txt file.
    # df_simplex = pd.DataFrame({'Simplex': simplex, 'Log-likelihood': simplex_llh, 'Frequency': simplex_frequency,
    #                            'Count': simplex_count, 'Average volume': simplex_average_v, 'Average tetrahedrality':
    #                                simplex_average_t})
    # df_aa = pd.DataFrame({'aa': aa, 'Frequency': aa_frequency, 'Count': aa_count})
    # with pd.ExcelWriter(os.path.join(llh_dir, 'simplex_llh_cutoff_12.xlsx')) as writer:
    #     df_simplex.to_excel(writer, sheet_name='Simplex')
    #     df_aa.to_excel(writer, sheet_name='aa')
    # writer.save()

    g = open(os.path.join(llh_dir, 'simplex_aa_llh_cutoff12.txt'), 'a')
    for i in range(len(simplex)):
        g.writelines(simplex[i] + '\t' + str("%.7f" % simplex_llh[i]) + '\t' + str("%.7f" % simplex_frequency[i]) +
                     '\t' + str(simplex_count[i]) + '\t' + str("%.3f" % simplex_average_v[i]) + '\t' +
                     str("%.3f" % simplex_average_t[i]) + '\n')
    g.close()

    g = open(os.path.join(llh_dir, 'aa_frequency_cutoff12.txt'), 'a')
    for i in range(len(aa)):
        g.writelines(aa[i] + '\t' + str("%.4f" % aa_frequency[i]) + '\t' + str(aa_count[i]) + '\n')
    g.close()


def edge_permutation_factor(edge):
    """Calculate permutation factor of each edge."""
    if edge[0:2] == edge[2:4]:
        return 1
    else:
        return 2


def calculate_edge_log_likelihood(tes_path, tes_files, analysis_path):
    """Calculate tetrahedron edges log-likelihood."""
    codon_dic = {}
    edge_dic = {}
    edge_length_dic = {}

    for tes_file in tes_files:
        f_tes = list(open(os.path.join(tes_path, "%s_tessellation" % tes_file), 'r'))

        res_pair_list = []
        for i in range(len(f_tes)):
            simplex = f_tes[i][:8]
            res = []
            for j in range(4):
                res.append([m.split('\t', j+3)[j+2] for m in [f_tes[i][:-1]]][0])

            len_list = []
            for j in range(6):
                len_list.append(float([m.split('\t', j+8)[j+7] for m in [f_tes[i][:-1]]][0]))

            for j in range(3):
                for k in range(j + 1, 4):
                    # Remove repeat edges in same file.
                    if [res[j], res[k]] not in res_pair_list:
                        # Add 1 to value if codon or edge appears.
                        check_key_in_dic(codon_dic, [simplex[j * 2:j * 2 + 2], simplex[k * 2:k * 2 + 2]])
                        # Sorted edge, make sure BA = AB.
                        sorted_edge_list = sorted([simplex[j * 2:j * 2 + 2], simplex[k * 2:k * 2 + 2]])
                        sorted_edge = str(sorted_edge_list[0] + sorted_edge_list[1])
                        check_key_in_dic(edge_dic, [sorted_edge])

                        # Add edge length to edge_length_dic.
                        if sorted_edge not in edge_length_dic:
                            if j == 0:
                                edge_length_dic[sorted_edge] = [len_list[k - 1]]
                            elif j == 1:
                                edge_length_dic[sorted_edge] = [len_list[k + 1]]
                            else:
                                edge_length_dic[sorted_edge] = [len_list[5]]
                        else:
                            if j == 0:
                                edge_length_dic[sorted_edge].append(len_list[k - 1])
                            elif j == 1:
                                edge_length_dic[sorted_edge].append(len_list[k + 1])
                            else:
                                edge_length_dic[sorted_edge].append(len_list[5])

                        res_pair_list.append([res[j], res[k]])

    # Calculate the frequency of each codon, edge.
    codon_frequency_dic = {k: v / total for total in (sum(codon_dic.values()),) for k, v in codon_dic.items()}
    edge_frequency_dic = {k: v / total for total in (sum(edge_dic.values()),) for k, v in edge_dic.items()}

    # Calculate log-likelihood for each edge.
    edge_llh_dic = {}
    for key in edge_frequency_dic:
        edge_llh_dic[key] = math.log10(edge_frequency_dic[key] / (edge_permutation_factor(key) * codon_frequency_dic[
            key[:2]] * codon_frequency_dic[key[2:4]]))

    # Calculate the average of edge length.
    edge_average_length_dic = {}
    for key, values in edge_length_dic.items():
        edge_average_length_dic[key] = statistics.mean(values)

    # Write into xlsx file.
    # df_edge_llh = pd.DataFrame(edge_llh_dic.items(), columns=['edge', 'Log-likelihood'])
    # df_codon_frequency = pd.DataFrame(codon_frequency_dic.items(), columns=['Codon', 'Frequency'])
    # df_edge_frequency = pd.DataFrame(edge_frequency_dic.items(), columns=['edge', 'Frequency'])
    # df_edge_count = pd.DataFrame(edge_dic.items(), columns=['edge', 'Count'])
    # df_codon_count = pd.DataFrame(codon_dic.items(), columns=['Codon', 'Count'])
    #
    # with pd.ExcelWriter(os.path.join(analysis_path, 'edge_codon_frequency_log_likelihood.xlsx')) as writer:
    #     df_edge_llh.to_excel(writer, sheet_name='edge log-likelihood')
    #     df_codon_frequency.to_excel(writer, sheet_name='Codon frequency')
    #     df_edge_frequency.to_excel(writer, sheet_name='edge frequency')
    #     df_edge_count.to_excel(writer, sheet_name='edge count')
    #     df_codon_count.to_excel(writer, sheet_name='Codon count')
    # writer.save()

    # Write into xlsx file with Edges and Codon.
    edge = []
    edge_llh = []
    edge_frequency = []
    edge_count = []
    edge_average_length = []
    codon = []
    codon_count = []
    codon_frequency = []
    for key in edge_llh_dic:
        edge.append(key)
        edge_llh.append(edge_llh_dic[key])
        edge_frequency.append(edge_frequency_dic[key])
        edge_count.append(edge_dic[key])
        edge_average_length.append(edge_average_length_dic[key])
    for key in codon_frequency_dic:
        codon.append(key)
        codon_count.append(codon_dic[key])
        codon_frequency.append(codon_frequency_dic[key])

    df_edge = pd.DataFrame({'Edges': edge, 'Log-likelihood': edge_llh, 'Frequency': edge_frequency, 'Count': edge_count,
                            'Average length': edge_average_length})
    df_codon = pd.DataFrame({'Codon': codon, 'Frequency': codon_frequency, 'Count': codon_count})
    with pd.ExcelWriter(os.path.join(analysis_path, 'short100_protein_edge_codon_llh.xlsx')) as writer:
        df_edge.to_excel(writer, sheet_name='Edge')
        df_codon.to_excel(writer, sheet_name='Codon')
    writer.save()


def identical_pair_counting(llh_dir, llh_file):
    f = pd.read_excel(os.path.join(llh_dir, llh_file), sheet_name='Edge')

    edges = f['Edges']
    rank = f['Rank']
    identical = []

    for i in range(len(edges)):
        if edges[i][:2] == edges[i][2:4]:
            identical.append(int(rank[i]))

    lim = [10, 50, 100, 200, 500]
    for i in range(len(lim)):
        count = 0
        for j in range(len(identical)):
            if identical[j] <= lim[i]:
                count += 1

        print("Top %s has: %s identical pairs." % (str(lim[i]), str(count)))


if __name__ == '__main__':
    simplex_analysis_dir = 'INDIVIDUAL/simplex_analysis'
    log_likelihood_dir = 'INDIVIDUAL/Log_likelihood'
    pdb_codon_dssp_dir = 'INDIVIDUAL/PDB_codon_dssp'
    tes_dir = 'INDIVIDUAL/PDB_tessellation'

    tes_files = []
    for filename in glob.glob(os.path.join(tes_dir, "*_tessellation")):
        tes_files.append(filename[73:])

    for cut_off in [12, 8, 4.8]:
        write_summary_simplex(tes_dir, tes_files, log_likelihood_dir, cut_off)

    sum_simplex_file = 'sum_all_tessellation_cutoff_8'
    calculate_simplex_log_likelihood(log_likelihood_dir, sum_simplex_file)
    calculate_simplex_log_likelihood_aa(log_likelihood_dir, sum_simplex_file)

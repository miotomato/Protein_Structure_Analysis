# Python 3.7

"""
File:           sql_db.py
Author:         Shengyuan Wang
Date:           Jun 1, 2020
Last update:    Jun 4, 2020

Purpose:        Create individual and protein-protein interaction sql databases, including description of each
                individual protein chain, protein sequence, nucleotide sequence, secondary structure, protein-protein
                interaction structure, bond length and type.

Input files:    Protein individual chain information, protein_pdb_codon_dssp files, ppi_codon files.
Output files:   Individual chain description, Individual chain structure, PPI chain description, PPI chain structure,
                PPI bond structure.
"""


import os
import glob
import pandas as pd
from math import isnan
import sqlite3
from sqlite3 import Error


def create_connection(path):
    connection = None
    try:
        connection = sqlite3.connect(path)
        print("Connection to SQLite DB successful")
    except Error as e:
        print(f"The error '{e}' occurred")

    return connection


def execute_query(connection, query):
    cursor = connection.cursor()
    try:
        cursor.execute(query)
        connection.commit()
        print("Query executed successfully")
    except Error as e:
        print(f"The error '{e}' occurred")


def codon_to_nt(codon):
    dic = {'A1': 'GCT', 'A2': 'GCC', 'A3': 'GCA', 'A4': 'GCG', 'C1': 'TGT', 'C2': 'TGC', 'D1': 'GAT', 'D2': 'GAC',
           'E1': 'GAA', 'E2': 'GAG', 'F1': 'TTT', 'F2': 'TTC', 'G1': 'GGT', 'G2': 'GGC', 'G3': 'GGA', 'G4': 'GGG',
           'H1': 'CAT', 'H2': 'CAC', 'I1': 'ATT', 'I2': 'ATC', 'I3': 'ATA', 'K1': 'AAA', 'K2': 'AAG', 'L1': 'TTA',
           'L2': 'TTG', 'L3': 'CTT', 'L4': 'CTC', 'L5': 'CTA', 'L6': 'CTG', 'M1': 'ATG', 'N1': 'AAT', 'N2': 'AAC',
           'P1': 'CCT', 'P2': 'CCC', 'P3': 'CCA', 'P4': 'CCG', 'Q1': 'CAA', 'Q2': 'CAG', 'R1': 'CGT', 'R2': 'CGC',
           'R3': 'CGA', 'R4': 'CGG', 'R5': 'AGA', 'R6': 'AGG', 'S1': 'TCT', 'S2': 'TCC', 'S3': 'TCA', 'S4': 'TCG',
           'S5': 'AGT', 'S6': 'AGC', 'T1': 'ACT', 'T2': 'ACC', 'T3': 'ACA', 'T4': 'ACG', 'V1': 'GTT', 'V2': 'GTC',
           'V3': 'GTA', 'V4': 'GTG', 'W1': 'TGG', 'Y1': 'TAT', 'Y2': 'TAC'}

    return dic[codon]


def read_ori_csv(ori_dir, ori_file):
    ori_dic = {'ID': [], 'loc': []}
    f_ori = pd.read_csv(os.path.join(ori_dir, ori_file), header=0)

    for i in range(len(f_ori['PDB ID'])):
        ori_dic['ID'].append(f_ori['PDB ID'][i] + f_ori['Chain ID'][i])
        ori_dic['loc'].append(i)

    return ori_dic


def pdb_codon_dssp_file(pdb_codon_dssp_dir, file, ori_dir, ori_file, ori_file_dic):
    """
    collect individual chain description and structure information
    """
    f_ori = pd.read_csv(os.path.join(ori_dir, ori_file), header=0)

    file_info = {}

    file_info['id'] = file[-20:-16]
    file_info['chain'] = file[-16]

    aa_seq = ''
    nt_seq = ''
    ss_seq = ''

    f = list(open(os.path.join(pdb_codon_dssp_dir, file), 'r'))

    start = int([m.split('\t', 3)[2] for m in [f[0]]][0])
    end = int([m.split('\t', 3)[2] for m in [f[-1]]][0])
    file_info['residue_num'] = str(start) + '-' + str(end)
    file_info['length'] = end - start + 1

    for line in f:
        aa_seq += line[0]
        nt_seq += codon_to_nt(line[:2])
        ss_seq += line[-2]

    file_info['aa_seq'] = aa_seq
    file_info['nt_seq'] = nt_seq
    file_info['ss_seq'] = ss_seq

    for i in range(len(ori_file_dic['ID'])):
        if ori_file_dic['ID'][i] == file[-20:-15]:
            file_info['resolution'] = f_ori['Resolution'][i]
            file_info['r_free'] = f_ori['R Free'][i]
            file_info['str_title'] = f_ori['Structure Title'][i]
            file_info['exp_method'] = f_ori['Exp. Method'][i]
            file_info['classification'] = f_ori['Classification'][i]
            file_info['uniprot'] = f_ori['Uniprot Acc'][i]
            file_info['source'] = f_ori['Source'][i]
            file_info['ec'] = f_ori['EC No'][i]

    for key in file_info:
        try:
            if isnan(file_info[key]):
                file_info[key] = 'NA'
        except TypeError:
            pass

    return file_info


def sql_insert_description(db_dir, db, file_info):
    """
    create individual chain description database
    """
    con = sqlite3.connect(os.path.join(db_dir, db))
    cur = con.cursor()

    try:
        cur.execute("insert into Individual_Chain_Description values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (file_info['id'], file_info['chain'], file_info['resolution'], file_info['r_free'],
                     file_info['str_title'], file_info['exp_method'], file_info['classification'],
                     file_info['uniprot'], file_info['source'], file_info['ec'], file_info['residue_num'],
                     file_info['length'], file_info['aa_seq'], file_info['nt_seq'], file_info['ss_seq']))
        con.commit()
        print("Query executed successfully")
    except Error as e:
        print(f"The error '{e}' occurred")


def sql_insert_structure(con, input_dir, input_file):
    """
    create individual chain structure database
    """
    cur = con.cursor()

    table_name = input_file[:5]
    query = "create table if not exists \'" + table_name + '\'' + \
            '(Codon TEXT, Amino_Acid TEXT, Residue_Num INT, cor_X FLOAT, cor_Y FLOAT, cor_Z FLOAT, ' \
            'Secondary_Structure Text);'
    cur.execute(query)

    f = list(open(os.path.join(input_dir, input_file), 'r'))

    for line in f:
        codon = line[:2]
        aa = [m.split('\t', 2)[1] for m in [line]][0]
        res_num = int([m.split('\t', 3)[2] for m in [line]][0])
        cor_x = float([m.split('\t', 4)[3] for m in [line]][0])
        cor_y = float([m.split('\t', 5)[4] for m in [line]][0])
        cor_z = float([m.split('\t', 6)[5] for m in [line]][0])
        ss = line[-2]

        try:
            query = "insert into \'" + table_name + '\'' + "values (?, ?, ?, ?, ?, ?, ?)"
            cur.execute(query, (codon, aa, res_num, cor_x, cor_y, cor_z, ss))
            con.commit()
        except Error as e:
            print("Error in:", input_file)
            print(f"The error '{e}' occurred")


def insert_sql_ppi(con, input_dir, input_file):
    """
    create protein-protein interface database
    """
    cur = con.cursor()

    table_name = input_file[:8]
    query = "create table if not exists \'" + table_name + '\'' + \
            '(Codon1 TEXT, Residue_Num1 INT, Secondary_Structure1 TEXT, Atom1 Text, cor1_X TEXT, cor1_Y TEXT, ' \
            'cor1_Z TEXT, Chain1 TEXT, Codon2 TEXT, Residue_Num2 INT, Secondary_Structure2 TEXT, Atom2 Text, ' \
            'cor2_X TEXT, cor2_Y TEXT, cor2_Z TEXT, Chain2 TEXT, Distance FLOAT, Bond TEXT);'
    cur.execute(query)

    f = list(open(os.path.join(input_dir, input_file), 'r'))

    for line in f:
        c1 = line[:2]
        res_num1 = int([m.split('\t', 2)[1] for m in [line]][0])
        ss1 = [m.split('\t', 3)[2] for m in [line]][0]
        atom1 = [m.split('\t', 4)[3] for m in [line]][0]
        cor1_x = [m.split('\t', 5)[4] for m in [line]][0]
        cor1_y = [m.split('\t', 6)[5] for m in [line]][0]
        cor1_z = [m.split('\t', 7)[6] for m in [line]][0]
        chain1 = [m.split('\t', 8)[7] for m in [line]][0]

        c2 = [m.split('\t', 10)[9] for m in [line]][0]
        res_num2 = int([m.split('\t', 11)[10] for m in [line]][0])
        ss2 = [m.split('\t', 12)[11] for m in [line]][0]
        atom2 = [m.split('\t', 13)[12] for m in [line]][0]
        cor2_x = [m.split('\t', 14)[13] for m in [line]][0]
        cor2_y = [m.split('\t', 15)[14] for m in [line]][0]
        cor2_z = [m.split('\t', 16)[15] for m in [line]][0]
        chain2 = [m.split('\t', 17)[16] for m in [line]][0]

        dis = float([m.split('\t', 18)[17] for m in [line]][0])
        bond = line[-2]

        try:
            query = "insert into \'" + table_name + '\'' + "values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, " \
                                                           "?, ?)"
            cur.execute(query, (c1, res_num1, ss1, atom1, cor1_x, cor1_y, cor1_z, chain1, c2, res_num2, ss2, atom2,
                                cor2_x, cor2_y, cor2_z, chain2, dis, bond))
            con.commit()
        except Error as e:
            print("Error in:", input_file)
            print(f"The error '{e}' occurred")


def insert_sql_tessellation(con, input_dir, input_file):
    """
    create protein tessellation database
    """
    cur = con.cursor()

    table_name = input_file[:5]
    query = "create table if not exists \'" + table_name + '\'' + \
            '(Simplex TEXT, ordered_Simplex TEXT, Residue_num1 INT, Residue_num2 INT, Residue_num3 INT, ' \
            'Residue_num4 INT, Residue_Distance INT, Edge1 FLOAT, Edge2 FLOAT, Edge3 FLOAT, Edge4 FLOAT, Edge5 FLOAT, ' \
            'Edge6 FLOAT, Volume FLOAT, Tetrahedrality FLOAT, SS1 TEXT, SS2 TEXT, SS3 TEXT, SS4 TEXT);'
    cur.execute(query)

    f = list(open(os.path.join(input_dir, input_file), 'r'))

    for line in f:
        sx = line[:8]
        ordered_sx = [m.split('\t', 2)[1] for m in [line]][0]

        res1 = int([m.split('\t', 3)[2] for m in [line]][0])
        res2 = int([m.split('\t', 4)[3] for m in [line]][0])
        res3 = int([m.split('\t', 5)[4] for m in [line]][0])
        res4 = int([m.split('\t', 6)[5] for m in [line]][0])
        res_dis = int([m.split('\t', 7)[6] for m in [line]][0])

        edge1 = float([m.split('\t', 8)[7] for m in [line]][0])
        edge2 = float([m.split('\t', 9)[8] for m in [line]][0])
        edge3 = float([m.split('\t', 10)[9] for m in [line]][0])
        edge4 = float([m.split('\t', 11)[10] for m in [line]][0])
        edge5 = float([m.split('\t', 12)[11] for m in [line]][0])
        edge6 = float([m.split('\t', 13)[12] for m in [line]][0])

        vol = float([m.split('\t', 14)[13] for m in [line]][0])
        th = float([m.split('\t', 15)[14] for m in [line]][0])

        ss1 = [m.split('\t', 16)[15] for m in [line]][0]
        ss2 = [m.split('\t', 17)[16] for m in [line]][0]
        ss3 = [m.split('\t', 18)[17] for m in [line]][0]
        ss4 = line[-2]

        try:
            query = "insert into \'" + table_name + '\'' + "values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, " \
                                                           "?, ?, ?)"
            cur.execute(query, (sx, ordered_sx, res1, res2, res3, res4, res_dis, edge1, edge2, edge3, edge4, edge5,
                                edge6, vol, th, ss1, ss2, ss3, ss4))
            con.commit()
        except Error as e:
            print("Error in:", input_file)
            print(f"The error '{e}' occurred")


def extract_sql_data(con, table_name, col_name):
    """
    extract data from database
    """
    con.row_factory = lambda cursor, row: row[0]
    cur = connection.cursor()

    x_sql = []
    query = "SELECT %s FROM %s" % (col_name, str('\'' + table_name + '\''))
    x_sql.append(cur.execute(query).fetchall())

    return x_sql


if __name__ == '__main__':
    pdb_codon_dssp_dir = "INDIVIDUAL/PDB_codon_dssp"
    pdb_ori_dir = "INDIVIDUAL/PDB_ori"
    ind_tes_dir = "INDIVIDUAL/PDB_tessellation"
    ppi_ori_dir = "PP/PDB_ori"
    ppi_codon_dir = "PP/PPI_codon_final"

    # create protein individual chain databases and tables
    connection = create_connection(os.path.join(pdb_ori_dir, 'PDB.db'))
    create_description_table = """
    CREATE TABLE IF NOT EXISTS Individual_Chain_Description (
        PDB_ID TEXT,
        Chain_ID TEXT,
        Resolution FLOAT,
        R_Free FLOAT,
        Structure_Title TEXT,
        Experiment_Method TEXT,
        Classification TEXT,
        Uniprot_Acc TEXT,
        Source TEXT,
        EC_No TEXT,
        PDB_Residue TEXT,
        AA_Chain_Length INT,
        AA_Seq TEXT,
        Nt_Seq TEXT,
        Secondary_Structure TEXT);
    """
    execute_query(connection, create_description_table)

    ori_file_dic = read_ori_csv(pdb_ori_dir, 'PDB_pc30_R2.0_Rf0.25_100_1000.csv')

    pdb_codon_dssp_files = []
    for filename in glob.glob(os.path.join(pdb_codon_dssp_dir, "*_pdb_codon_dssp")):
        pdb_codon_dssp_files.append(filename[-20:])

    for codon_file in pdb_codon_dssp_files:
        pdb_file_info = pdb_codon_dssp_file(pdb_codon_dssp_dir, codon_file, pdb_ori_dir,
                                            'PDB_pc30_R2.0_Rf0.25_100_1000.csv', ori_file_dic)
        sql_insert_description(pdb_ori_dir, 'PDB.db', pdb_file_info)

    connection = create_connection(os.path.join(pdb_ori_dir, 'Ind_structure.db'))
    pdb_codon_dssp_files = []
    for filename in glob.glob(os.path.join(pdb_codon_dssp_dir, "*_pdb_codon_dssp")):
        pdb_codon_dssp_files.append(filename[-20:])

    for codon_file in pdb_codon_dssp_files:
        sql_insert_structure(connection, pdb_codon_dssp_dir, codon_file)

    # create protein-protein interfaces database and tables
    connection = create_connection(os.path.join(ppi_ori_dir, 'PPI_structure.db'))
    ppi_codon_files = []
    for filename in glob.glob(os.path.join(ppi_codon_dir, "*_PPI_codon_final")):
        ppi_codon_files.append(filename[-24:])

    for ppi_file in ppi_codon_files:
        insert_sql_ppi(connection, ppi_codon_dir, ppi_file)

    # create tessellation database
    connection = create_connection(os.path.join(pdb_ori_dir, 'Ind_tessellation.db'))
    tes_files = []
    for filename in glob.glob(os.path.join(ind_tes_dir, "*_tessellation")):
        tes_files.append(filename[-18:])

    for tes_file in tes_files:
        insert_sql_tessellation(connection, ind_tes_dir, tes_file)


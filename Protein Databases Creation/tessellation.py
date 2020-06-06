# Python 3.7

"""
File:           Delaunay_tessellation
Author:         Shengyuan Wang
Date:           Jun 5, 2019
Last update:    Jun 4, 2020

Purpose:        Find protein residue neighbors using Delaunay Tessellation.

Input files:    Protein pdb_codon_dssp files.
Output files:   Protein tessellation files. Tessellation analysis.

Notices:        1. All Delaunay simplex whose edges greater than 12 Angstrom will be removed.
                2. Output Delaunay tessellation following characters: Simplex, sorted_simplex, residue numbers,
                residue diatance, edge lengths, volume, tetrahedrality, DSSP.
                3. Simplices are classified based on counting of consecutive number of residue numbers.
                4. Density plots are performed for volumes and tetrahedarlity for each class of simplices.
"""


import os
import glob
import math
from pyhull.delaunay import DelaunayTri
import matplotlib.pyplot as plt
import seaborn as sns


def distance(a, b):
    """Calculate the edge length of tetrahedron."""
    dist = math.sqrt((float(a[0]) - float(b[0])) ** 2 + (float(a[1]) - float(b[1])) ** 2 +
                     (float(a[2]) - float(b[2])) ** 2)
    return str('%.3f' % dist)


def determinant_3x3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
            m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1]) +
            m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]))


def subtract(a, b):
    return (a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2])


def tetrahedron_calc_volume(a, b, c, d):
    """Calculate the volume of tetrahedron."""
    return str('%.3f' % (abs(determinant_3x3((subtract(a, b), subtract(b, c), subtract(c, d),))) / 6.0))


def tetrahedron_calc_tetrahedrality(l1, l2, l3, l4, l5, l6):
    """Calculate the tetrahedrality of tetrahedron."""
    list = [l1,l2,l3,l4,l5,l6]
    list = [float(i) for i in list]
    average = sum(list) / len(list)
    sum_edge = 0
    for i in range(0,5):
        for j in range(i+1,6):
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

    return l1, l2, l3, l4, l5, l6, volume, tetrahedrality



def tessellation_file(input_dir, input_file, cut_distance, output_dir):
    """
    Write tessellation file based on pdb_codon_dssp files, all simplex whose edges greater than distance_cut will be
    removed.
    """
    f = list(open(os.path.join(input_dir, '%s_pdb_codon_dssp' % input_file), 'r'))
    g = open(os.path.join(output_dir, "%s_tessellation" % input_file[:5]), 'a')

    codon_list = []
    residue_num = []
    coordinates = []
    dssp_list = []

    # Extract codon, residue number, coordinates, dssp info.
    for i in range(len(f)):
        codon = f[i][:2]
        res = [m.split('\t', 3)[2] for m in [f[i][:-1]]][0]
        cor_x = float([m.split('\t', 4)[3] for m in [f[i][:-1]]][0])
        cor_y = float([m.split('\t', 5)[4] for m in [f[i][:-1]]][0])
        cor_z = float([m.split('\t', 6)[5] for m in [f[i][:-1]]][0])
        dssp = [m.split('\t', 6)[6] for m in [f[i][:-1]]][0]

        codon_list.append(codon)
        coordinates.append([cor_x, cor_y, cor_z])
        residue_num.append(res)
        dssp_list.append(dssp)

    # Perform Delaunay Tessellation using pyhull, inputs are coordinates.
    tri = DelaunayTri(coordinates)

    # Sort DelaunayTri.vertices inside list first, then entire list.
    sort_tri_vertices = []
    for i in range(len(tri.vertices)):
        cache0_tri = sorted(tri.vertices[i])
        sort_tri_vertices.append(cache0_tri)
    sort_tri_vertices = sorted(sort_tri_vertices)

    for i in range(len(sort_tri_vertices)):
        cache_tri = sort_tri_vertices[i]

        # Four neighbor codons as a simplex.
        simplex = codon_list[cache_tri[0]]+codon_list[cache_tri[1]]+codon_list[cache_tri[2]]+codon_list[cache_tri[3]]

        # Rewrite simplex based on alphabet.
        simplex_sort = sorted([codon_list[cache_tri[0]], codon_list[cache_tri[1]], codon_list[cache_tri[2]],
                               codon_list[cache_tri[3]]])
        sorted_simplex = simplex_sort[0]+simplex_sort[1]+simplex_sort[2]+simplex_sort[3]

        # Write residue number based on simplex (not sorted simplex).
        r1 = residue_num[cache_tri[0]]
        r2 = residue_num[cache_tri[1]]
        r3 = residue_num[cache_tri[2]]
        r4 = residue_num[cache_tri[3]]

        # Calculate residue distance.
        residue_distance = int(r4)-int(r1)-3
        cor1 = coordinates[cache_tri[0]]
        cor2 = coordinates[cache_tri[1]]
        cor3 = coordinates[cache_tri[2]]
        cor4 = coordinates[cache_tri[3]]

        # Calculate tetrahedron edges length, volume and tetrahedrality.
        l1, l2, l3, l4, l5, l6, volume, tetrahedrality = tetrahedron_calculation(cor1, cor2, cor3, cor4)

        # Write dssp based on simplex (not sorted simplex).
        ss1 = dssp_list[cache_tri[0]]
        ss2 = dssp_list[cache_tri[1]]
        ss3 = dssp_list[cache_tri[2]]
        ss4 = dssp_list[cache_tri[3]]

        edge_list = [l1, l2, l3, l4, l5, l6]
        edge_list = [float(j) for j in edge_list]

        # Check if all edge length less or equal to 12 Angstrom.
        # Write simplex, sorted_simplex, residue numbers, residue distance, edges length, volume, tetrahedrality, dssp.
        if all(edge <= cut_distance for edge in edge_list):
            g.writelines(simplex + '\t' + sorted_simplex + '\t' + r1 + '\t' + r2 + '\t' + r3 + '\t' + r4 + '\t' +
                         str(residue_distance) + '\t' + l1 + '\t' + l2 + '\t' + l3 + '\t' + l4 + '\t' + l5 + '\t' + l6 +
                         '\t' + volume + '\t' + tetrahedrality + '\t' + ss1 + '\t' + ss2 + '\t' + ss3 + '\t' + ss4 +
                         '\n')

    g.close()


def check_consecutive_number(res1, res2, res3, res4):
    # 4 consecutive residues.
    if res4 - res3 == res3 - res2 == res2 - res1 == 1:
        return 4
    # 3 consecutive residues.
    elif res3 - res2 == res2 - res1 == 1 and res4 - res3 != 1:
        return 3
    elif res4 - res3 == res3 - res2 == 1 and res2 - res1 != 1:
        return 3
    # 2 pairs 2 consecutive residues.
    elif res2 - res1 == res4 - res3 == 1 and res3 - res2 != 1:
        return 2
    # All separate residues.
    elif res2 - res1 != 1 and res3 - res2 != 1 and res4 - res3 != 1:
        return 0
    # 1 pair 2 consecutive residues.
    else:
        return 1


def classify_simplex(input_dir, tes_file, output_dir):
    """Classify simplices based on the number of consecutive residues"""
    f = list(open(os.path.join(input_dir, tes_file), 'r'))
    g_0 = open(os.path.join(output_dir, 'sum_simplex_class1111'), 'a')     # All separate residues.
    g_1 = open(os.path.join(output_dir, 'sum_simplex_class211'), 'a')      # 1 pair 2 consecutive residues.
    g_2 = open(os.path.join(output_dir, 'sum_simplex_class22'), 'a')       # 2 pairs 2 consecutive residues.
    g_3 = open(os.path.join(output_dir, 'sum_simplex_class31'), 'a')       # 3 consecutive residues.
    g_4 = open(os.path.join(output_dir, 'sum_simplex_class4'), 'a')        # 4 consecutive residues.

    for i in range(len(f)):
        simplex = f[i][:8]
        simplex_list = [simplex[:2], simplex[2:4], simplex[4:6], simplex[6:8]]
        sorted_simplex = [m.split('\t', 2)[1] for m in [f[i][:-1]]][0]
        res1 = int([m.split('\t', 3)[2] for m in [f[i][:-1]]][0])
        res2 = int([m.split('\t', 4)[3] for m in [f[i][:-1]]][0])
        res3 = int([m.split('\t', 5)[4] for m in [f[i][:-1]]][0])
        res4 = int([m.split('\t', 6)[5] for m in [f[i][:-1]]][0])
        res_list = [res1, res2, res3, res4]
        res_dis = int([m.split('\t', 7)[6] for m in [f[i][:-1]]][0])
        volume = [m.split('\t', 14)[13] for m in [f[i][:-1]]][0]
        tetrahedrality = [m.split('\t', 15)[14] for m in [f[i][:-1]]][0]

        # Simplify dssp into 1 alphabet, N = 'NA', M = ' ', G = 3-turn helix, H = 4-turn helix, I = 5-turn helix,
        # T = hydrogen bond turn, E = beta-sheet, B = residue in isolated beta-bridge, S = bend, C = coil.
        dssp_list = []
        for j in range(4):
            dssp_list.append([m.split('\t', 16+j)[15+j] for m in [f[i][:-1]]][0])

        # Sorted res_list and dssp_clean_list based on simplex _list.
        sorted_simplex_res = sorted(zip(simplex_list, res_list))
        sorted_simplex_dssp = sorted(zip(simplex_list, dssp_list))

        # Classify simplex.
        check_result = check_consecutive_number(res1, res2, res3, res4)
        simplex_class = [g_0, g_1, g_2, g_3, g_4]

        # Write into different simplex class files.
        simplex_class[check_result].writelines(tes_file[:5] + '\t' + sorted_simplex + '\t' +
                                               str(sorted_simplex_res[0][1]) + '\t' + str(sorted_simplex_res[1][1]) +
                                               '\t' + str(sorted_simplex_res[2][1]) + '\t' +
                                               str(sorted_simplex_res[3][1]) + '\t' + str(res_dis) + '\t' + volume +
                                               '\t' + tetrahedrality + '\t' + sorted_simplex_dssp[0][1] + '\t'
                                               + sorted_simplex_dssp[1][1] + '\t' + sorted_simplex_dssp[2][1] + '\t' +
                                               sorted_simplex_dssp[3][1] + '\n')


def density_plot(file_path, input_files, input, xlim_set):
    """Draw density plot for input = volume/tetrahedrality. Need input xlim for plot."""
    files = []
    lists = [[], [], [], [], []]
    labels = []
    for file in input_files:
        files.append(list(open(os.path.join(file_path, file), 'r')))
        labels.append(file[12:])

    loc = int
    if input == 'volume':
        loc = 7
    elif input == 'tetrahedrality':
        loc = 8

    for i in range(len(files)):
        for j in range(len(files[i])):
            need_plot = float([m.split('\t', loc+1)[loc] for m in [files[i][j][:-1]]][0])
            lists[i].append(need_plot)

    for i in range(len(lists)):
        # kde, i.e. Kernel density estimation is more computationally, but it makes the plot more smooth.
        sns.distplot(lists[i], hist=False, kde=True, label=labels[i])

    plt.xlabel(input)
    plt.ylabel('Density')
    plt.xlim(left=0, right=xlim_set)
    # After plt.show() is called, a new figure is created. To deal with this, save the figure first then display.
    fig = plt.gcf()
    plt.show()
    fig.savefig(os.path.join(file_path, 'Density plot of %s of simplices in 5 different classes.eps' % input))
    plt.close(fig)


if __name__ == '__main__':
    # Calculate Delaunay tessellation for all pdb_codon_dssp files.
    pdb_codon_dssp_dir = 'INDIVIDUAL/PDB_codon_dssp'
    tes_dir = 'INDIVIDUAL/PDB_tessellation'
    cut_distance = 12
    simplex_analysis_dir = 'INDIVIDUAL/simplex_analysis'

    for filename in glob.glob(os.path.join(pdb_codon_dssp_dir, "*_pdb_codon_dssp")):
        tessellation_file(pdb_codon_dssp_dir, filename[-20:-15], cut_distance, tes_dir)

    # Classify simplices based on counts of consecutive residues.
    for filename in glob.glob(os.path.join(tes_dir, "*_tessellation")):
        classify_simplex(tes_dir, filename[-18:], simplex_analysis_dir)

    # Calculate percentage of secondary structure in different simplex classes.
    sum_simplex = []
    for filename in glob.glob(os.path.join(simplex_analysis_dir, "sum_simplex_class*")):
        sum_simplex.append(filename[73:])

    sum_simplex = sorted(sum_simplex, reverse=True)
    density_plot(simplex_analysis_dir, sum_simplex, 'volume', 80)
    density_plot(simplex_analysis_dir, sum_simplex, 'tetrahedrality', 0.6)

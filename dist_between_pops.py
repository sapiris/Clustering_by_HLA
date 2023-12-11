import sys
import os
import numpy as np
from collections import Counter
import pickle
import shutil
import pathlib
import pandas as pd
import sys



sys.path.insert(0, os.path.join(".."))
import math
import gzip
# from validation.calc_KL import calc_KL as KL

import json
import argparse


def calc_AF_dist(file_real_freq, em_file):
    dict_haps = {}
    dict_alleles = {"A": {}, "B": {}, "C": {}, "DRB1": {}, "DQB1": {}}

    with open(file_real_freq) as true_file:
        # with gzip.open(file_real_freq, 'rb') as zf:
        # true_file = [x.decode('utf8').strip() for x in zf.readlines()]
        for line in true_file:
            hap, pop, freq = line.strip().split(',')
            hap = ('~').join(sorted(hap.split('~')))
            dict_haps[hap] = [float(freq), 0]
            for allele in hap.split('~'):
                locus = allele.split('*')[0]
                if allele not in dict_alleles[locus]:
                    dict_alleles[locus][allele] = [0, 0]
                dict_alleles[locus][allele][0] += float(freq)

    # with gzip.open(em_file, 'rb') as zf:
    with open(em_file) as true_file:
        # true_file = [x.decode('utf8').strip() for x in zf.readlines()]
        for line in true_file:
            line = line.strip().split(',')
            if line[0] in dict_haps:
                dict_haps[line[0]][1] = (float(line[2]))
            else:
                dict_haps[line[0]] = [0, float(line[2])]
            hap, freq = line[0], line[2]
            for allele in hap.split('~'):
                locus = allele.split('*')[0]
                if allele not in dict_alleles[locus]:
                    dict_alleles[locus][allele] = [0, 0]
                dict_alleles[locus][allele][1] += float(freq)

    sum = 0
    for key in dict_haps:
        p = dict_haps[key][0]
        q = dict_haps[key][1]
        sum += math.pow((p - q), 2)

    list_d = []
    all_d = 0
    for locus in dict_alleles:
        d = 0
        for key in dict_alleles[locus]:
            d += (dict_alleles[locus][key][0] * dict_alleles[locus][key][1]) ** 0.5
        d = 2 / math.pi * math.sqrt(2 * (1 - d))
        list_d.append(d)
        all_d = + d ** 2

    D = math.sqrt(all_d)

    return D

def calc_euclidain_dist_by_haplo(file_real_freq, em_file):
    dict_haps = {}
    dict_alleles = {}

    with open(file_real_freq) as true_file:
        # with gzip.open(file_real_freq, 'rb') as zf:
        # true_file = [x.decode('utf8').strip() for x in zf.readlines()]
        for line in true_file:
            hap, pop, freq = line.strip().split(',')
            hap = ('~').join(sorted(hap.split('~')))
            dict_haps[hap] = [float(freq), 0]
            """for allele in hap.split('~'):
                if allele not in dict_alleles:
                    dict_alleles[allele] = [0, 0]
                dict_alleles[allele][0] += float(freq)"""

    # with gzip.open(em_file, 'rb') as zf:
    with open(em_file) as true_file:
        # true_file = [x.decode('utf8').strip() for x in zf.readlines()]
        for line in true_file:
            line = line.strip().split(',')
            if line[0] in dict_haps:
                dict_haps[line[0]][1] = (float(line[2]))
            else:
                dict_haps[line[0]] = [0, float(line[2])]
            """hap, freq = line[0], line[2]
            for allele in hap.split('~'):
                if allele not in dict_alleles:
                    dict_alleles[allele] = [0, 0]
                dict_alleles[allele][1] += float(freq)"""

    sum = 0
    for key in dict_haps:
        p = dict_haps[key][0]
        q = dict_haps[key][1]
        sum += math.pow((p - q), 2)

    """sum_allele = 0
    for key in dict_alleles:
        p = dict_alleles[key][0]
        q = dict_alleles[key][1]
        sum_allele += math.pow((p - q), 2)"""
    return  math.sqrt(sum)#math.sqrt(sum_allele)  # ,

def calc_euclidain_dist(file_real_freq, em_file):
    dict_haps = {}
    dict_alleles = {}

    with open(file_real_freq) as true_file:
        # with gzip.open(file_real_freq, 'rb') as zf:
        # true_file = [x.decode('utf8').strip() for x in zf.readlines()]
        for line in true_file:
            hap, pop, freq = line.strip().split(',')
            hap = ('~').join(sorted(hap.split('~')))
            dict_haps[hap] = [float(freq), 0]
            for allele in hap.split('~'):
                if allele not in dict_alleles:
                    dict_alleles[allele] = [0, 0]
                dict_alleles[allele][0] += float(freq)

    # with gzip.open(em_file, 'rb') as zf:
    with open(em_file) as true_file:
        # true_file = [x.decode('utf8').strip() for x in zf.readlines()]
        for line in true_file:
            line = line.strip().split(',')
            if line[0] in dict_haps:
                dict_haps[line[0]][1] = (float(line[2]))
            else:
                dict_haps[line[0]] = [0, float(line[2])]
            hap, freq = line[0], line[2]
            for allele in hap.split('~'):
                if allele not in dict_alleles:
                    dict_alleles[allele] = [0, 0]
                dict_alleles[allele][1] += float(freq)

    """sum = 0
    for key in dict_haps:
        p = dict_haps[key][0]
        q = dict_haps[key][1]
        sum += math.pow((p - q), 2)"""

    sum_allele = 0
    for key in dict_alleles:
        p = dict_alleles[key][0]
        q = dict_alleles[key][1]
        sum_allele += math.pow((p - q), 2)
    return math.sqrt(sum_allele)  # , math.sqrt(sum)


def calc_dist(pops_list, simulation, path_freqs, path_output, byalleles = True, ed_dist = True):
    pops_len = len(pops_list)
    len_list = []

    for race in pops_list:
       len_list.append(sum(1 for line in open('output/freqs/freq_' + race + '.csv')))

    set_cala = set()
    dist_matrix = np.zeros(([pops_len, pops_len]))

    for i in range(pops_len):
        # print(i)
        race1 = pops_list[i]
        f1 = f"{path_freqs}/hpf_allele_{race1}.csv"
        # f1 = f'../imputation/graph_generation/data/2014/{race1}.freqs.gz'
        for j in range(i + 1, pops_len):
            race2 = pops_list[j]
            print(f"{race1},{race2}")
            set_cala.add(f"{race1},{race2}")
            set_cala.add(f"{race2},{race1}")
            f2 = f"{path_freqs}/hpf_allele_{race2}.csv"
            # f2 =  f'../imputation/graph_generation/data/2014/{race2}.freqs.gz'
            if byalleles and ed_dist:
                euc_dist = calc_euclidain_dist(f1, f2)
            """elif ed_dist and not byalleles:
                euc_dist = calc_euclidain_dist_by_haplo(f1, f2)
            elif  not byalleles and not ed_dist:
                kl1 = calc_KL(len_list[i], f1, f2)
                kl2 = calc_KL(len_list[j], f2, f1)
                euc_dist = (kl1 + kl2) / 2
            elif  byalleles and not ed_dist:
                kl1 = calc_KL_alleles(len_list[i], f1, f2)
                kl2 = calc_KL_alleles(len_list[j], f2, f1)
                euc_dist = (kl1 + kl2) / 2"""

            dist_matrix[i, j] = euc_dist
            dist_matrix[j, i] = euc_dist

    d = pd.DataFrame(dist_matrix, index=pops_list, columns=pops_list)

    d.to_excel(f'{path_output}/dist_matrix/dist_matrix_{simulation}.xlsx')  ##check end

    pickle.dump(dist_matrix, open(f'{path_output}/dist_matrix/dist_matrix_{simulation}.pkl', "wb"))


    #return d



def main(simulation, path_to_freqs, f_races_size_path, path_output):
    list_include = []
    with open(f_races_size_path) as size_f:
        for line in size_f:
            line = line.strip().split(',')
            #if float(line[1]) > 50:
            list_include.append(line[0])


    #cut_off = 50  # clusters pops just with number of individuals bigger than this cut_off
    print(simulation)

    calc_dist(list_include, simulation, path_to_freqs, path_output, byalleles = True, ed_dist = True)







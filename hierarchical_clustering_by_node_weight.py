import pandas as pd
import pickle
import copy
import networkx as nx
import math
import numpy as np

def clusters_size(f_races_size_path):
    dict_pop_size = {}
    sum_all =0
    #with open("../clustering/race_Size.csv") as size_file:
    #with open("../EM/data/Israel/race_count.csv") as size_file:
    #with open("../EM/data/indai_race_size.csv") as size_file:
    with open(f_races_size_path) as size_file:
        for line in size_file:
            line = line.strip().split(',')
            dict_pop_size[line[0]] =  float(line[1])
            sum_all += float(line[1])

    dict_pop_size = { key:val/sum_all for key, val in dict_pop_size.items()}

    return dict_pop_size




def ward(dist_matrix, pops_list, num_clusters, dict_pop_size):
    d = pd.DataFrame(dist_matrix, index=pops_list, columns=pops_list)
    d = d.replace(0, 99)
    clusters = list(d.columns)
    while len(clusters) > num_clusters:
        dict_dist = {}

        min_value = d.min().min()  # Get the minimum value in the entire DataFrame
        min_value_indices = d[d == min_value].stack().index  # Get the indices of the minimum value

        # Extract the column and row names
        race1 = min_value_indices[0][0]
        race2 = min_value_indices[0][1]
        min_pair = [race1, race2]

        list_dist = []
        c_i = min_pair[0]
        c_j = min_pair[1]
        n_i =  1#dict_pop_size[c_i]
        n_j = 1#dict_pop_size[c_j]
        for c_k in clusters:
            if c_k not in min_pair:
                n_k = 1# dict_pop_size[c_k]
                dist = (n_i + n_k)/(n_i + n_j + n_k) * d.loc[c_i, c_k] + \
                       (n_j + n_k) / (n_i + n_j + n_k) * d.loc[c_j, c_k] - \
                       (n_k)/(n_i + n_j + n_k) * d.loc[c_i, c_j]

                list_dist.append(dist)#dict_dist[race1 + '+' + race2]


        d = d.drop(labels=min_pair,  axis=0)
        d = d.drop(labels=min_pair, axis=1)
        dict_pop_size[f"{min_pair[0]}-{min_pair[1]}"] =  dict_pop_size[min_pair[1]] + dict_pop_size[min_pair[0]]
        del dict_pop_size[min_pair[0]]
        del dict_pop_size[min_pair[1]]
        d[f"{min_pair[0]}-{min_pair[1]}"] = list_dist
        list_dist.append(99)
        d.loc[f"{min_pair[0]}-{min_pair[1]}"] = list_dist
        clusters = list(d.columns)

        #for race in d.columns:
        #    for r in race1.split("-"):
    for c in clusters:
        print(c)



def find_closest_cluster(d, dict_pop_size, list_pop_to_add, graph_clusters, clusters):
    #dict_dist = {}
    min_dist = 1000000000000000
    min_pair = []

    for pop in list_pop_to_add:
        for race1 in clusters:
            if race1 != pop:
                if d.loc[race1, pop] < min_dist:
                    min_dist =  d.loc[race1, pop]
                    min_pair = [race1, pop]

    print(min_pair[1])
    if "Brazil" == min_pair[1]:
         k = 0
    list_pop_to_add.remove(min_pair[1])
    list_dist = []
    c_i = min_pair[0]
    c_j = min_pair[1]
    n_i = dict_pop_size[c_i]
    n_j = dict_pop_size[c_j]
    for c_k in d.columns:
        if c_k in clusters and c_k not in min_pair:
            n_k = dict_pop_size[c_k]
            dist = (n_i + n_k) / (n_i + n_j + n_k) * d.loc[c_i, c_k] + \
                   (n_j + n_k) / (n_i + n_j + n_k) * d.loc[c_j, c_k] - \
                   (n_k) / (n_i + n_j + n_k) * d.loc[c_i, c_j]

            list_dist.append(dist)  # dict_dist[race1 + '+' + race2]
        else:
            list_dist.append(99)

    #d = d.drop(labels=min_pair, axis=0)
    #d = d.drop(labels=min_pair, axis=1)
    new_cluster = f"{min_pair[0]}-{min_pair[1]}"
    dict_pop_size[new_cluster] = dict_pop_size[min_pair[1]] + dict_pop_size[min_pair[0]]
    graph_clusters.add_node(new_cluster, size=dict_pop_size[new_cluster], dist = min_dist)
    graph_clusters.add_edge(new_cluster, min_pair[0])
    graph_clusters.add_edge(new_cluster, min_pair[1])
    del dict_pop_size[min_pair[0]]
    del dict_pop_size[min_pair[1]]
    if f"{min_pair[0]}-{min_pair[1]}" not in d:
        d[f"{min_pair[0]}-{min_pair[1]}"] = list_dist
        list_dist.append(99)
        d.loc[f"{min_pair[0]}-{min_pair[1]}"] = list_dist
    clusters.remove(min_pair[0])
    clusters.append(new_cluster)

    return clusters



def ward_with_min_size(dist_matrix, pops_list, num_clusters, dict_pop_size, f_cluster_path, threshold_for_cluster = 1/10):# 120000):
    threshold_for_cluster = 1/num_clusters * threshold_for_cluster
    graph_clusters = nx.DiGraph()
    for pop,size in dict_pop_size.items():
        graph_clusters.add_node(pop, size = size)
    list_d = []
    d = pd.DataFrame(dist_matrix, index=pops_list, columns=pops_list)
    list_d.append(d.copy())
    d = d.replace(0, 99)
    d_all = copy.deepcopy(d)
    clusters = list(d.columns)
    dict_pop_size_all =copy.deepcopy(dict_pop_size)
    while len(clusters) > num_clusters:
        dict_dist = {}
        min_dist = 1000000000000000
        min_pair = []

        min_value = d.min().min()  # Get the minimum value in the entire DataFrame
        min_value_indices = d[d == min_value].stack().index  # Get the indices of the minimum value

        # Extract the column and row names
        race1 = min_value_indices[0][0]
        race2 = min_value_indices[0][1]
        min_pair = [race1, race2]

        list_dist = []
        c_i = min_pair[0]
        c_j = min_pair[1]
        n_i =  dict_pop_size[c_i]
        n_j = dict_pop_size[c_j]
        for c_k in clusters:
            if c_k not in min_pair:
                n_k = dict_pop_size[c_k]
                dist = (n_i + n_k)/(n_i + n_j + n_k) * d.loc[c_i, c_k] + \
                       (n_j + n_k) / (n_i + n_j + n_k) * d.loc[c_j, c_k] - \
                       (n_k)/(n_i + n_j + n_k) * d.loc[c_i, c_j]

                list_dist.append(dist)#dict_dist[race1 + '+' + race2]


        d = d.drop(labels=min_pair,  axis=0)
        d = d.drop(labels=min_pair, axis=1)
        new_cluster = f"{min_pair[0]}-{min_pair[1]}"
        dict_pop_size[new_cluster] =  dict_pop_size[min_pair[1]] + dict_pop_size[min_pair[0]]
        graph_clusters.add_node(new_cluster, size=dict_pop_size[new_cluster], dist= min_value)
        graph_clusters.add_edge(new_cluster,min_pair[0] )
        graph_clusters.add_edge(new_cluster, min_pair[1])
        del dict_pop_size[min_pair[0]]
        del dict_pop_size[min_pair[1]]
        d[f"{min_pair[0]}-{min_pair[1]}"] = list_dist
        list_dist.append(99)
        d.loc[f"{min_pair[0]}-{min_pair[1]}"] = list_dist
        clusters = list(d.columns)
        list_d.append([d.copy(),copy.deepcopy(graph_clusters), copy.deepcopy(dict_pop_size)])

        dict_pop_size_all[new_cluster] =  dict_pop_size_all[min_pair[1]] + dict_pop_size_all[min_pair[0]]
        list_dist = []
        c_i = min_pair[0]
        c_j = min_pair[1]
        n_i =  dict_pop_size_all[c_i]
        n_j = dict_pop_size_all[c_j]
        for c_k in list(d_all.columns):
            if c_k not in min_pair:
                n_k = dict_pop_size_all[c_k]
                dist = (n_i + n_k)/(n_i + n_j + n_k) * d_all.loc[c_i, c_k] + \
                       (n_j + n_k) / (n_i + n_j + n_k) * d_all.loc[c_j, c_k] - \
                       (n_k)/(n_i + n_j + n_k) * d_all.loc[c_i, c_j]

                list_dist.append(dist)#dict_dist[race1 + '+' + race2]
            else:
                list_dist.append(99)
        d_all[f"{min_pair[0]}-{min_pair[1]}"] = list_dist
        list_dist.append(99)
        d_all.loc[f"{min_pair[0]}-{min_pair[1]}"] = list_dist

        #for race in d.columns:
        #    for r in race1.split("-"):
    #for c in clusters:
    #    print(c)
    for cluster_round in range(len(list_d)-1,0,-1):
        d_idx, graph_clusters_idx, dict_pop_size_idx = list_d[cluster_round]
        clusters, list_pop_to_add = [], []
        for pop, size in dict_pop_size_idx.items():
            if size < threshold_for_cluster:
                list_pop_to_add.append(pop)
            else:
                clusters.append(pop)
        if len(clusters) == num_clusters:
            final_clusters = clusters
            while list_pop_to_add:
                final_clusters = find_closest_cluster(d_all, dict_pop_size_idx, list_pop_to_add, graph_clusters_idx, final_clusters)

            """print("\n\n new*******************")

            for c in final_clusters:
                print(c)"""

            f_o = open(f_cluster_path, "w")
            f_o.write("Cluster,Races\n")
            for i, c in enumerate(final_clusters):
                f_o.write(f"{i},{c.replace('-', ';')}\n")
                print(c)

            return final_clusters


def main(num_clusters, f_races_size_path, f_cluster_path, simulation, path_output):
    dist_matrix = pickle.load(
        open(f"{path_output}/dist_matrix/dist_matrix_{simulation}.pkl", "rb"))  # f'output/pkl/dist_matrix_AF_{simulation}.pkl'
    dict_pop_size = clusters_size(f_races_size_path)

    ward_with_min_size(dist_matrix, list(dict_pop_size.keys()), num_clusters, dict_pop_size, f_cluster_path, 1 / 10)







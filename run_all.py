from freq_by_typing import estimate_freqs
import hierarchical_clustering_by_node_weight
from step_3_for_small_races import assign_small_races
import os
import pathlib
import dist_between_pops
import json
from EM.run_em import run_em_def

#step 1 - create clusters
conf_file = "conf/minimal-em-configuration.json"
with open(conf_file) as f:
    json_conf = json.load(f)

typing_path = json_conf.get("input_file")

simulation = json_conf.get("dataset_name")
output_path = json_conf.get("clustering_output_path")
num_clusters =  json_conf.get("num_clusters") #10,20,40
f_cluster_path = f"{output_path}/{num_clusters}clusters_{simulation}.csv"
#if dir not exist create it
if not os.path.exists(output_path):
    os.makedirs(output_path)

path_to_freqs = f"{output_path}/freqs_{simulation}"
pathlib.Path(path_to_freqs).mkdir(parents=False, exist_ok=True)

f_races_size_path = f"{output_path}/{simulation}_race_size_over_threshold.csv"

estimate_freqs(typing_path, f_races_size_path, path_to_freqs,  size_threshold = 5)
dist_between_pops.main(simulation, path_to_freqs, f_races_size_path, output_path)

hierarchical_clustering_by_node_weight.main(num_clusters, f_races_size_path, f_cluster_path, simulation, output_path)




#change races in input files to clusters
def change_races_to_clusters(f_cluster_path, typing_path, output_typing):
    f_update_file = open(output_typing, "w")
    race_cluster_dict = {}
    with open(f_cluster_path, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            cluster, races = line.strip().split(",")
            for race in races.split(";"):
                race_cluster_dict[race] = f"Cluster_{cluster}"

    with open(typing_path, "r") as f:
        for line in f.readlines():
            line = line.strip().split(",")
            for idx in [2,3]:
                list_races = []
                for race in line[idx].split(";"):
                    if race in race_cluster_dict:
                        list_races.append(race_cluster_dict[race])
                    else:
                        list_races.append(race)
                line[idx] = ";".join(list_races)
            f_update_file.write(",".join(line) + "\n")
    f_update_file.close()

change_races_to_clusters(f_cluster_path, typing_path, output_typing=json_conf.get("imputation_in_file"))

#step 2 - run EM
run_em_def(conf_file)

#step 3 - assign small races to clusters
assign_small_races(conf_file, num_clusters, output_path, simulation)

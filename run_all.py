from freq_by_typing import estimate_freqs
import hierarchical_clustering_by_node_weight
import os
import pathlib
import dist_between_pops

typing_path = ""
simulation = "US"
output_path = "output/"
#if dir not exist create it
if not os.path.exists(output_path):
    os.makedirs(output_path)

path_to_freqs = f"{output_path}/freqs_{simulation}"
pathlib.Path(path_to_freqs).mkdir(parents=False, exist_ok=True)

f_races_size_path = f"{output_path}/{simulation}_race_size_over_threshold.csv"

estimate_freqs(typing_path, f_races_size_path, path_to_freqs,  size_threshold = 50)
dist_between_pops.main(simulation, path_to_freqs, f_races_size_path, output_path)

for num_clusters in [5]: #10,20,40
    f_cluster_path = f"output/{num_clusters}clusters_{simulation}.csv"
    hierarchical_clustering_by_node_weight.main(num_clusters, f_races_size_path, f_cluster_path, simulation, output_path)

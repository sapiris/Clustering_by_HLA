from freq_by_typing import estimate_freqs
import hierarchical_clustering_by_node_weight
import os
import pathlib
import dist_between_pops
import argparse
from EM.run_em import run_em_def
from step_3_for_small_races import assign_small_races
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the analysis on a CSV file.")
    parser.add_argument("--input_file", type=str, help="Path to the typing input CSV", default="input/input_example.csv")
    parser.add_argument("--output", type=str, default="output_example", help="The output folder")
    parser.add_argument("--cluster_size_threshold", type=int, default=6)
    parser.add_argument("--min_pop_size", type=int, default=0)
    parser.add_argument("--num_clusters", type=int, default=2)
    parser.add_argument("--conf_file", type=str, default="conf/minimal-em-configuration.json")
    args = parser.parse_args()
    typing_path = args.input_file
    output_path = args.output
    #if dir not exist create it
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    path_to_freqs = f"{output_path}/freqs"
    pathlib.Path(path_to_freqs).mkdir(parents=False, exist_ok=True)

    f_races_size_path = f"{output_path}/race_size_over_threshold.csv"

    is_all_pop_over_threshold = estimate_freqs(typing_path, f_races_size_path, path_to_freqs,  size_threshold = args.min_pop_size)
    dist_between_pops.main( path_to_freqs, f_races_size_path, output_path)

    num_clusters = args.num_clusters
    f_cluster_path = f"{output_path}/{num_clusters}clusters.csv"
    hierarchical_clustering_by_node_weight.main(num_clusters, f_races_size_path, f_cluster_path, output_path, args.cluster_size_threshold)

    if not is_all_pop_over_threshold:

        # change races in input files to clusters
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
                    for idx in [2, 3]:
                        list_races = []
                        for race in line[idx].split(";"):
                            if race in race_cluster_dict:
                                list_races.append(race_cluster_dict[race])
                            else:
                                list_races.append(race)
                        line[idx] = ";".join(list_races)
                    f_update_file.write(",".join(line) + "\n")
            f_update_file.close()


        conf_file = args.conf_file
        with open(conf_file) as f:
            json_conf = json.load(f)
        change_races_to_clusters(f_cluster_path, typing_path, output_typing=json_conf.get("imputation_in_file"))

        # step 2 - run EM
        run_em_def(conf_file)

        # step 3 - assign small races to clusters
        assign_small_races(conf_file, num_clusters, output_path, typing_path)


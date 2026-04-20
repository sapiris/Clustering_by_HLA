import json

def assign_small_races(config_file, num_of_cluster, output_clustering, simulation):
        with open(config_file) as f:
            json_conf = json.load(f)

        output_dir = json_conf.get("imputation_out_path", "output") + "/"
        config = {
            "imputation_input_file": json_conf.get("input_file"),#imputation_in
            "imputation_out_hap_freq_file": output_dir + json_conf.get("imputation_out_hap_freq_filename"),
            "pops": json_conf.get("populations"),
            # "KL_freq_file_to_compare": json_conf.get("KL_file_to_compare"),
            # "KL_real_pop_size": json_conf.get("KL_real_size")
        }


        dict_occurance_of_race = {}

        dict_id_races = {}
        with open(config["imputation_input_file"]) as input_f:#ToDo add imputation in file
            for line in input_f:
                line = line.strip().split(',')
                races = line[2] + ';' + line[3]
                races = races.split(';')
                races[:] = (value for value in races if value != '')
                dict_id_races[line[0]] = races
                for race in races:
                    dict_occurance_of_race[race] = dict_occurance_of_race.get(race, 0) + 1/len(races)

        input_f.close()
        dict_all_orig_races = {}
        list_finel_races = config["pops"] #todo cahnge races list
        with open(config["imputation_out_hap_freq_file"]) as imputat_res:#todo change results file
            res = imputat_res.readline()
            while res:
                id = res.split(",")[0]
                sum = 0
                list_person = [0]*num_of_cluster

                # find the all results of single person
                # sum the prob of single person
                while (res and res.split(",")[0] == id):
                    res = res.split(",")
                    p = float(res[3])
                    sum += p
                    p = p/2
                    race1 = res[1].split(';')[1]
                    race2 = res[2].split(';')[1]
                    list_person[list_finel_races.index(race1)] += p
                    list_person[list_finel_races.index(race2)] += p

                    res = imputat_res.readline()
                #normalise individual prob sum to 1
                list_person = [value/sum for value in list_person]
                id_orig_races = dict_id_races[id]

                if len(id_orig_races) > 0:
                    #divide sum to number of individual races in SIRE(so total prob of individual will be 1)
                    list_person = [value / len(id_orig_races) for value in list_person]

                    for race in id_orig_races:
                        if not race in dict_all_orig_races:
                            dict_all_orig_races[race] = [0]*num_of_cluster
                        for i in range(num_of_cluster):

                            dict_all_orig_races[race][i] += (list_person[i])
        imputat_res.close()

        number_of_races = 0
        list_all_probs = [0]*num_of_cluster
        for race in dict_all_orig_races:
            race_size = dict_occurance_of_race[race]
            number_of_races += race_size
            for i in range(len(list_all_probs)):
                list_all_probs[i] = list_all_probs[i] + race_size*dict_all_orig_races[race][i]

        list_all_probs = [value/number_of_races for value in list_all_probs]
        for race in dict_all_orig_races:
            for i in range(len(dict_all_orig_races[race])):
                dict_all_orig_races[race][i] /= list_all_probs[i]


        dict_my_map = {}
        """for pop in list_finel_races:
            dict_my_map[pop] = pop"""
        with open(f'{output_clustering}/{num_of_cluster}clusters_{simulation}.csv') as race_cluster: #todo add name file
            for line in race_cluster:
                if "Cluster" in line:
                    continue
                line = line.strip().split(',')
                list_races = line[1].split(';')
                for race in list_races:
                    dict_my_map[race] = 'Cluster_' + str(line[0])
        race_cluster.close()


        f_res = open(f'{output_clustering}/check_pop_clusters_IL_{simulation}cluster_ward.csv', 'w')
        f_res.write('race_by_EzMi,')
        for i in range(num_of_cluster):
            f_res.write(f'prob_to_Cluster_{i},')
        f_res.write ('cluster with max prob, my tag,match\n')
        count = 0
        for key,value in dict_all_orig_races.items():
            line = key + ','
            for p in value:
                line += str(p)
                line += ','
            idx = value.index(max(value))
            line += list_finel_races[idx]
            line += ','
            if key in dict_my_map:
                line += dict_my_map[key]

                if dict_my_map[key] == list_finel_races[idx]:
                    line += ',1'
                    count+=1
                else:
                    line += ',0'
            else:
                line += ','
            f_res.write(line + '\n')
        f_res.close()
        print(count)
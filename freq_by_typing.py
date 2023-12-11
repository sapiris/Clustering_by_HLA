import sys
import os

sys.path.insert(0, os.path.join(".."))
#from utils import open_gl_string, clean_up_gl
from grim.imputation.impute import Imputation, clean_up_gl
def add_haps_to_dict(dict_loci, phases_list, num_races):
    sum = 0
    freq_allele = 0
    for phases in phases_list:
        for j in phases:
            for i in j:
                for k in i:
                    freq_allele += 1
    freq_allele *= num_races
    for phases in phases_list:
        for j in phases:
            for i in j:
                for k in i:
                    phase = ('~').join(k)
                    if not phase in dict_loci:
                        # if len(dict_loci) > self.max_size:
                        # dict_loci = cut_dict(dict_loci, self.min_size)
                        dict_loci[phase] = (2 / freq_allele)
                    else:
                        dict_loci[phase] += (2 / freq_allele)
                    sum += (2 / freq_allele)
    return sum  # , dict_loci


def write_freqs_to_file(dict_loci, sum, pop, freq_file):
    f_freq = open(freq_file, 'w')
    for hap in dict_loci:
        f_freq.write(hap + ',' + pop + ',' + str(dict_loci[hap] / sum) + '\n')

    f_freq.close()


def is_full_haplo(hap):
    if all(loci in hap for loci in ["A", "B", "C", "DRB1", "DQB1"]):
        return True
    return False


def estimate_freqs(typing_path, f_races_size_path, path_to_freqs,  size_threshold = 50):
    haps_count = {}
    imputation = Imputation()
    dict_race_loci = {}
    #read file to count the pop sizes
    dict_race_size = {}

    with open(typing_path) as lines:
        for name_gl in lines:
            name_gl = name_gl.strip().split(',')
            races = name_gl[2] + ';' + name_gl[3]
            race_list = races.split(';')
            for race in race_list:
                dict_race_size[race] = dict_race_size.get(race, 0) + 1/len(race_list)

    f_races_size = open(f_races_size_path, "w")
    list_include = []
    for r, size in dict_race_size.items():
        if size > size_threshold:
            dict_race_loci[r] = {}
            haps_count[r] = 0
            list_include.append(r)
            f_races_size.write(f"{r},{size}\n")

    num_samples = 0
    dict_race_count = {}

    with open(typing_path) as lines:
        for name_gl in lines:
            if ',' in name_gl:
                name_gl = name_gl.strip().split(',')
            else:
                name_gl = name_gl.strip().split('%')
            name_gl[1] = clean_up_gl(name_gl[1])
            # if name_gl[2] != name_gl[3]:
            #    continue
            race = name_gl[2] + ';' + name_gl[3]
            race_list = race.split(';')
            if is_full_haplo(name_gl[1]):

                phases_list = imputation.open_gl_string(name_gl[1], 1000)
                #phases_list = open_gl_string(name_gl[1], 1000)
                if not phases_list:
                    continue

                for race in race_list:
                    if race in list_include:
                        dict_race_count[race] = dict_race_count.get(race, 0) + 1 / len(race_list)
                        haps_count[race] += add_haps_to_dict(dict_race_loci[race], phases_list, len(race_list))

                num_samples += 1
    print("num_samples_from_freqs", num_samples)
    lines.close()

    """f_oo = open("data/Israel/race_size_by_typing2.csv", "w")
    for race, size in dict_race_count.items():
        f_oo.write(f"{race},{size}\n")"""
    # calculate the probability of each haplotype and write to file
    for r, dict_loci in dict_race_loci.items():
        write_freqs_to_file(dict_loci, haps_count[r], r,f"{path_to_freqs}/hpf_allele_{r}.csv")



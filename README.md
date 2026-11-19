# HAPLOCLUST

HAPLOCLUST automatically clusters self-defined groups using HLA typing data to minimize false disease associations and maximize within-group genetic similarity.

## Input

Provide a CSV file containing HLA typing data for all populations. Each row represents one individual and contains an ID, HLA typing in GL String format, and two self-identified race and ethnicity fields (`SIRE1` and `SIRE2`).

See [`input/input_example.csv`](input/input_example.csv) for an example.

## Run HAPLOCLUST

Run the commands below from the repository root. Replace the example input path and output directory with your own.

### With uv

With [uv](https://docs.astral.sh/uv/getting-started/installation/) installed, run:

```bash
uv run --locked python run_all.py --input_file input/input_example.csv --output results
```

`uv` sets up the project environment and installs the dependencies recorded in `uv.lock`.

### With Python and pip

In an activated Python environment, install the dependencies and run the script:

```bash
python -m pip install -r requirements.txt
python run_all.py --input_file input/input_example.csv --output results
```

## Parameters

| Parameter | Description |
| --- | --- |
| `--input_file` | Path to the input CSV containing HLA typing data. |
| `--output` | Path to the output directory. |
| `--cluster_size_threshold` | User-defined threshold used to calculate the minimum cluster size. |
| `--min_pop_size` | Minimum number of individuals a population needs to participate in the initial clustering. Smaller populations are assigned to clusters in a subsequent step. |
| `--num_clusters` | Requested number of clusters. |
| `--conf_file` | Path to the GRIMME configuration file used to assign populations smaller than `--min_pop_size`. |

## GRIMME configuration for small populations

The pipeline uses the GRIMME configuration only when at least one population has fewer individuals than `--min_pop_size`. In the configuration file:

- Include the cluster labels `Cluster_0` through `Cluster_<num_clusters - 1>` in the population list. For example, if `--num_clusters` is `3`, include `Cluster_0`, `Cluster_1`, and `Cluster_2`.
- Set `imputation_in_file` to the path where HAPLOCLUST will generate the input file for GRIMME. Use a different path from `--input_file` to avoid overwriting your original input CSV.
- Set `run_just_SR_EM` to `false`.

For the full list of configuration options, see the [GRIMME configuration documentation](https://github.com/nmdp-bioinformatics/py-graph-imputation/blob/master/conf/README.md).

## Output

HAPLOCLUST writes `<X>clusters.csv` to the directory specified by `--output`, where `<X>` is the value of `--num_clusters`.

Each row represents one cluster. The cluster number comes first, followed by a comma and the population names separated by semicolons.

If the small-population assignment step runs, HAPLOCLUST also writes `<X>clusters_witsmallpop.csv` to the same directory. This file includes all population assignments.

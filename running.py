import os

import drugstone
import time
import pandas as pd
import random

drugstone.set_api('https://api.drugst.one/')

arguments = {
    "drug": ["trustrank", "closeness", "degree", "proximity"],
    "drug-target": ["multisteiner"
        , "keypathwayminer", "trustrank", "closeness", "degree", "betweenness"
                    ]
}

test_files = {
    "10": [
        "./seed_files/amyloidosis-10seeds.tsv",
        "./seed_files/delirium-10seeds.tsv",
        "./seed_files/metabolic_disease-10seeds.tsv",
        "./seed_files/pseudohypoaldosteronism-10seeds.tsv",
        "./seed_files/sickle_cell_anemia-10seeds.tsv",
    ], "25": [
        "./seed_files/astrocytoma_(excluding_glioblastoma)-25seeds.tsv",
        "./seed_files/pancytopenia-26seeds.tsv",
        "./seed_files/colitis-25seeds.tsv",
        "./seed_files/thrombocytopenia-29seeds.tsv",
        "./seed_files/retinal_disorder-27seeds.tsv",
    ],
    "50": [
        "./seed_files/leukemia-55seeds.tsv",
        "./seed_files/medulloblastoma-50seeds.tsv",
        "./seed_files/psoriasis-57seeds.tsv",
        "./seed_files/nervous_system_disorder-53seeds.tsv",
        "./seed_files/ulcerative_colitis-63seeds.tsv"
    ],
    "100": [
        "./seed_files/HIV_infectious_disease-103seeds.tsv",
        "./seed_files/visual_epilepsy-101seeds.tsv",
        "./seed_files/retinitis_pigmentosa-104seeds.tsv",
        "./seed_files/cirrhosis_of_liver-103seeds.tsv",
        "./seed_files/cardiomyopathy-130seeds.tsv",
    ]
}



random.seed(42)
def read_tests(test_files):
    gene_lists = list()
    for case, file_list in test_files.items():
        target_size = int(case)
        for file in file_list:
            with open(file, 'r') as f:
                genes = list()
                for gene in f.readlines():
                    if gene.startswith("#"):
                        continue
                    genes.append(gene.strip())
                while len(genes) > target_size:
                    genes.remove(genes[random.randrange(len(genes) - 1)])
                gene_lists.append(genes)
    return gene_lists


seed_lists = read_tests(test_files)
parameters = {
    "identifier": "entrez",
    "algorithm": "trustrank",
    "target": "drug-target",
    "ppiDataset": 'NeDRex',
    "pdiDataset": "NeDRex",
}
drugstone.print_license()
drugstone.accept_license()


def run_task_and_get_runningtime(genes, parameters):
    start_time = time.time()
    try:
        task = drugstone.new_task(genes, parameters)
    except:
        print(f"Error with {genes} and {parameters}")
    end_time = time.time()
    return end_time - start_time


stats = list()

total_runs = (0 + len(arguments["drug-target"])) * len(seed_lists)
run_nr = 0

for k, v in arguments.items():
    for algorithm in v:
        for genes in seed_lists:
            run_nr += 1
            print(f"Run {run_nr} of {total_runs}")
            parameters["algorithm"] = algorithm
            parameters["target"] = k
            runtime = run_task_and_get_runningtime(genes, parameters)
            print(f"Runtime for {algorithm} is {runtime}")
            stats.append(
                {'algorithm': algorithm, 'target': k, 'runtime': runtime, 'seed_size': len(genes), 'seeds': genes, })
df = pd.DataFrame.from_records(stats)

df.to_csv('./drugstone_runtime.tsv', index=False, sep='\t')

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv("./drugstone_runtime.tsv", sep="\t")

sns.set(style="whitegrid")

plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 12
})

targets = ['drug', 'drug-target']
custom_ticks = [10, 25, 50, 100]

titles = {
    'drug': "Runtime Scaling for Drug Ranking",
    'drug-target': "Runtime Scaling for Drug-Target Identification"
}

algorithm_label_map = {
    "betweenness": "Betweenness centrality",
    "closeness": "Harmonic centrality",
    "degree": "Degree centrality",
    "keypathwayminer": "KeyPathwayMiner",
    "multisteiner": "MultiSteiner",
    "proximity": "Network proximity",
    "trustrank": "TrustRank"
}

data['algorithm'] = data['algorithm'].map(algorithm_label_map)

for target in targets:
    target_data = data[data['target'] == target]

    plt.figure(figsize=(10, 6))
    ax = plt.subplot(1, 1, 1)

    sns.lineplot(data=target_data, x='seed_size', y='runtime', hue='algorithm', style='algorithm',
                 markers=True, dashes=False, ci='sd', err_style='band', ax=ax)

    ax.set_ylabel('Runtime (seconds)', fontsize=14)
    ax.set_xlabel('Number of Seed Genes', fontsize=14)
    ax.set_title(titles[target], fontsize=16)
    ax.set_yscale('log')
    ax.set_xticks(custom_ticks)
    ax.set_xticklabels(custom_ticks)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend(title='Algorithm', title_fontsize='13', fontsize='11', loc='upper left')

    plt.tight_layout()
    plt.savefig(f'./drugstone_runtime_scaling_{target}.png', dpi=600)

    # plt.show()

# saved_plots = [f'./drugstone_runtime_scaling_{target}.png' for target in targets]
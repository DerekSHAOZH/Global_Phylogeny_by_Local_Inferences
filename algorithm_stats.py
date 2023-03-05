import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_algorithm_stats(filename):
    df = pd.read_csv (filename)
    df.plot(x="k", y=["MST nodes", "Global tree vertices", "V size","MP score of local", "Global tree edges"])
    plot = df.plot(x="k", y=["MST nodes", "Global tree vertices", "V size","MP score of local", "Global tree edges"])
    plt.xlabel("iteration (k)")
    filename_without_folder = filename.split('/')[-1]

    dataset_threshold = filename_without_folder.split('.')[0]
    dataset_threshold_parts = dataset_threshold.split('_')
    dataset = dataset_threshold_parts[0]
    threshold = dataset_threshold_parts[1]
    
    plt.title(f'influenza_{dataset}.fasta, threshold = {threshold}')
    fig = plot.get_figure()
    plot_filename = f'Result/{dataset_threshold}_stats_plot.png'
    print(f'plot_filename: {plot_filename}')
    fig.savefig(plot_filename)

if __name__ == "__main__":
    plot_algorithm_stats(sys.argv[1])
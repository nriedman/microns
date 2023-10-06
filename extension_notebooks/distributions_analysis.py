import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os


def plot_distribution(data, title, xlabel, ylabel, null=None, bins=50, plot_path=None):
    """
    Plot the distribution of the data.
    :param data: The data to plot.
    :param title: The title of the plot.
    :param xlabel: The label of the x axis.
    :param ylabel: The label of the y axis.
    :param bins: The number of bins to use.
    :return: None
    """
    fig, ax = plt.subplots()

    # Plot the histogram
    ax.hist(data, bins=bins, label='Data', density=True)

    if null:
        # Plot the null distribution
        ax.hist(null, bins=bins, label='Null', alpha=0.5, density=True)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Plot the mean, standard deviation, and standard error of the mean with legend
    mean = np.mean(data)
    std = np.std(data)
    sem = std / np.sqrt(len(data))

    ax.axvline(mean, color='r', linestyle='dashed', linewidth=1, label=f'Mean: {str(mean)[:5]} ± {str(sem)[:5]}')
    ax.axvline(mean + std, color='k', linestyle='dashed', linewidth=1, label=f"Standard Deviation: {str(std)[:5]}")
    ax.axvline(mean - std, color='k', linestyle='dashed', linewidth=1)

    ax.legend()

    if plot_path:
        output = os.path.join(plot_path, title.replace(' ', '_') + '_a_0')
        while os.path.exists(output + '.pkl'):
            output = output[:-1] + str(int(output[-1]) + 1)
        output += '.png'

        plt.savefig(output)
    else:
        plt.show()


def main(data_path, null_path=None, save_stats=False, plot_path=None):
    # Load the data
    # A dictionary keys for each sequence length and values are the list of lcs comparisons
    lcs_dict = pd.read_pickle(data_path)
    distributions = {}

    if null_path:
        null_dict = pd.read_pickle(null_path)
        null_distributions = {}

    # Preprocess the data
    for length, sequence_similarities in lcs_dict.items():
        print(f"Processing length {length}...")

        cur_distribution = distributions.setdefault(length, [])
        for seq, lcs_comp in sequence_similarities:
            cur_distribution.append(lcs_comp[0][1])
        
        print(f"Processed {len(cur_distribution)} sequences.", end='\n\n')
    
    if null_path:
        for length, sequence_similarities in null_dict.items():
            print(f"Processing null length {length}...")

            cur_distribution = null_distributions.setdefault(length, [])
            for seq, lcs_comp in sequence_similarities:
                cur_distribution.append(lcs_comp[0][1])
            
            print(f"Processed {len(cur_distribution)} null sequences.", end='\n\n')
    
    # Plot the distributions
    all_stats = {}
    for length, data in distributions.items():
        stats = all_stats.setdefault(length, {})
        stats['min'] = min(data)
        stats['max'] = max(data)
        stats['percentiles'] = np.percentile(data, [5, 10, 25, 50, 75, 90, 95])
        stats['mean'] = np.mean(data)
        stats['median'] = np.median(data)
        stats['std'] = np.std(data)
        stats['var'] = np.var(data)
        stats['num_seqs'] = len(data)

        bins = np.arange(min(data), max(data)+1, 1)
        print(f"Length {length}:")
        print(f"Min Similarity: {min(data)}")
        print(f"Max Similarity: {max(data)}")
        print(f"Mean Similarity: {np.mean(data)}")
        print(f"Median Similarity: {np.median(data)}")
        print(f"Standard Deviation: {np.std(data)}")
        print(f"Variance: {np.var(data)}")
        print(f"Number of Sequences: {len(data)}")
        print('_'*50)

        if null_path and length in null_distributions:
            plot_distribution(data, 'Similarity Distribution for short sequence length {} with Null'.format(length), 'Similarity', 'Frequency', null=null_distributions[length], bins=bins, plot_path=plot_path)
        else:
            plot_distribution(data, 'Similarity Distribution for short sequence length {}'.format(length), 'Similarity', 'Frequency', bins=bins, plot_path=plot_path)

    if save_stats:
        print("Saving stats to file...")

        # Get the path to the directory containing the data file
        pd.to_pickle(all_stats, data_path[:-4] + 'stats_.pkl')


if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser(description='Plot the distribution of the data.')
    parser.add_argument('--data', type=str, help='The data file to plot (pickle).')
    parser.add_argument('--null', type=str, default=None, help='The null data file to plot (pickle).')
    parser.add_argument('--save_stats', type=bool, default=False, help='Whether to save the stats to a file.')
    parser.add_argument('--plot_path', type=str, default=None, help='The path to save the plot to (png).')

    args = parser.parse_args()

    main(args.data, args.null, args.save_stats, args.plot_path)
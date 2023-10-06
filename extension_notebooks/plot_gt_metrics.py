
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import argparse
import utils as ut

def plot_metrics(metrics, title, save_path=None):
    # Given a dataframe of metrics, plot the mean precision and recall vs. synapse count with error bars
    # metrics: dataframe with columns ['synapse_count', 'precision', 'recall']
    # title: the title of the plot
    # save_path: the path to save the plot to
    # Return: None
    metrics = metrics.sort_values(by='synapse_count')
    synapse_counts = metrics['synapse_count'].unique()
    precision_means = []
    precision_stds = []
    precision_maxs = []
    precision_mins = []

    recall_means = []
    recall_stds = []
    recall_maxs = []
    recall_mins = []

    for synapse_count in synapse_counts:
        recall = metrics[metrics['synapse_count'] == synapse_count]['recall']
        precision = metrics[metrics['synapse_count'] == synapse_count]['precision']

        precision_means.append(float(precision.mean()))
        precision_stds.append(float(precision.std()))
        precision_maxs.append(float(precision.max()))
        precision_mins.append(float(precision.min()))

        recall_means.append(float(recall.mean()))
        recall_stds.append(float(recall.std()))
        recall_maxs.append(float(recall.max()))
        recall_mins.append(float(recall.min()))
    
    plt.errorbar(synapse_counts, precision_means, yerr=precision_stds, label='Precision')
    plt.errorbar(synapse_counts, recall_means, yerr=recall_stds, label='Recall')

    plt.fill_between(synapse_counts, precision_mins, precision_maxs, alpha=0.2, label='Precision Range')
    plt.fill_between(synapse_counts, recall_mins, recall_maxs, alpha=0.2, label="Recall Range")

    plt.xlabel('Synapse Count')
    plt.ylabel('Precision and Recall')
    plt.title(title)
    plt.legend()
    plt.show()

    if save_path != None:
        plt.savefig(save_path)
        print("Saved plot to {}".format(save_path))


def compile_metrics(dir):
    # Loop through all the directories in dir
    # For each directory, load the metrics and add them to a dataframe
    # Plot the dataframe
    metrics = pd.DataFrame(columns=['synapse_count', 'precision', 'recall'])
    print("Compiling metrics...")
    for i, subdir in enumerate(os.listdir(dir)):
        subdir_path = os.path.join(dir, subdir)
        if os.path.isdir(subdir_path):
            metrics_path = os.path.join(subdir_path, 'gt_metrics.pkl')
            metrics_dict = pd.read_pickle(metrics_path)

            # Find the synapse count in the subdir name, flagged by an s
            flags = subdir.split('_')
            synapse_count = int(flags[2][1:])

            data = {'synapse_count': [synapse_count],
                    'precision': [float(metrics_dict['precision'])],
                    'recall': [float(metrics_dict['recall'])]}
            
            metrics = pd.concat([metrics, pd.DataFrame(data=data)], ignore_index=True)
        
        print(f"Done {i+1}/{len(os.listdir(dir))}")

    print("Done compiling metrics")
    return metrics.apply(pd.to_numeric)


def main(dir, save_path, metrics_path=None):
    # Compile the metrics from all the samples in the dir
    if metrics_path == None:
        metrics = compile_metrics(dir)

        # Save the compiled metrics in the dir
        print("Saving compiled metrics...")
        metrics.to_csv(os.path.join(dir, 'compiled_metrics.csv'))
        print("Saved compiled metrics to {}".format(os.path.join(dir, 'compiled_metrics.csv')))
    else:
        metrics = pd.read_csv(metrics_path, index_col=0)
        print("Loaded metrics from {}".format(metrics_path))
    
    print("_"*80)

    # Plot the metrics
    print("Plotting metrics...")
    plot_metrics(metrics, 'Precision and Recall vs. Synapse Count', save_path=save_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot the ground truth metrics")
    parser.add_argument('--dir', type=str, help='The directory containing the ground truth samples')
    parser.add_argument('--metrics', type=str, default=None, help='The path to the compiled metrics')
    parser.add_argument('--save_path', type=str, default=None, help='The path to save the plot to')
    args = parser.parse_args()

    main(args.dir, args.save_path, metrics_path=args.metrics)
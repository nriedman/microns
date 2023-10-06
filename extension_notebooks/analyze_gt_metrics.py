
import numpy as np
import pandas as pd
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import utils as ut
import os

"""
Food for thought / future investigation:
    - Strange pattern observed in precision and recall: Both decrease from 0-10 synapse count, but precision increases 50-200 and recall continues to decrease
    - Happens with strict definition of true positive (exact match)
    - Come back and generate a distribution for different values of synapse_count
"""


"""
Want: Recall and Precision <==> False positives and False negatives

Recall: How many of the ground truth paths are in the exc paths
Precision: How many of the exc paths are in the ground truth paths

False positives: How many exc paths are not in the ground truth paths
False negatives: How many ground truth paths are not in the exc paths

Precision = 1 - false positives / total exc paths
Recall = 1 - false negatives / total ground truth paths

"""


def is_subsequence(subsequence, sequences):
    """
    Check if a subsequence is in a list of sequences
    """
    for sequence in sequences:
        if len(sequence) < len(subsequence):
            continue

        for i in range(len(sequence) - len(subsequence) + 1):
            if sequence[i:i+len(subsequence)] == subsequence:
                return True
    
    return False


def get_metrics(_dir, exc_synapse, cells_df, _ignore_burst=False):
    # Load relevant tables
    print("Loading {}".format(_dir + '/ground_truth_synapses.csv'))
    gt = pd.read_csv(_dir + '/ground_truth_synapses.csv', index_col=0)
    gt.set_index('synapse_id', inplace=True)

    print("Loading {}".format(_dir + '/participant_ids.pkl'))
    cell_ids = pd.read_pickle(_dir + '/participant_ids.pkl')

    print("Slimming down tables...")
    cells_df = cells_df.loc[cell_ids]
    exc_synapse = exc_synapse.loc[exc_synapse['post_pt_root_id'].isin(cell_ids), gt.columns]

    print("Loaded all tables")
    print("_"*80)

    # Generate msts for the exc and ground truth tables
    print("Generating MSTs...")
    exc_msts, _ = ut.generate_msts(exc_synapse, cells_df)
    gt_msts, _ = ut.generate_msts(gt, cells_df)

    # Extract the paths from the msts
    print("Extracting paths...")

    paths = {}
    for i, (gt_mst, exc_mst) in enumerate(zip(gt_msts, exc_msts)):
        cell_id_gt = gt_mst.graph['cell_id']
        cell_id_exc = exc_mst.graph['cell_id']
        assert cell_id_gt == cell_id_exc

        print(cell_id_exc, f'{i+1}/{len(gt_msts)}', end='\r')
        paths[cell_id_exc] = (ut.get_paths(gt_mst, root=-1, ignore_burst=_ignore_burst), ut.get_paths(exc_mst, root=-1, ignore_burst=_ignore_burst))
    
    print("Done extracting paths. Cleaning...")

    # Clean the paths
    paths_cleaned = {}
    for cell_id, (gt_paths, exc_paths) in paths.items():
        cleaned = paths_cleaned.setdefault(cell_id, ([], []))
        for gt_path, exc_path in zip(gt_paths, exc_paths):
            cleaned[0].append(gt_path) if len(gt_path) >= 3 else None
            cleaned[1].append(exc_path) if len(exc_path) >= 3 else None
    
    print("Done cleaning paths. Saving...")

    # Save the paths
    pd.to_pickle(paths_cleaned, _dir + '/gt_cleaned_paths.pkl')

    print("Saved to {}".format(_dir + '/gt_cleaned_paths.pkl'))
    print("_"*80)

    # Generate the metrics
    print("Generating metrics...")
    metrics = {}

    total_exc = 0
    total_gt = 0
    false_positives = 0
    false_negatives = 0

    for cell_id, (gt_paths, exc_paths) in paths_cleaned.items():
        total_exc += len(exc_paths)
        total_gt += len(gt_paths)

        for gt_path, exc_path in zip(gt_paths, exc_paths):
            false_positives += 0 if exc_path in gt_paths else 1
            false_negatives += 0 if is_subsequence(gt_path, exc_paths) else 1
    

    metrics['total_exc'] = total_exc
    metrics['total_gt'] = total_gt
    metrics['false_positives'] = false_positives
    metrics['false_negatives'] = false_negatives

    print(f'Total excitatory paths: {total_exc}')
    print(f'Total ground truth paths: {total_gt}')
    print(f'False positives: {false_positives}')
    print(f'False negatives: {false_negatives}')

    precision = 1 - false_positives / total_exc
    recall = 1 - false_negatives / total_gt

    print(f'Precision: {precision}')
    print(f'Recall: {recall}')

    metrics['precision'] = precision
    metrics['recall'] = recall


    print("Done generating metrics. Saving...")
    pd.to_pickle(metrics, _dir + '/gt_metrics.pkl')
    print("Saved to {}".format(_dir + '/gt_metrics.pkl'))


def main(_dir, exc_table='synapses_w_ids.csv', cell_table='cells_no_repeats.csv', _ignore_burst=False):
    # Load relevant tables

    print("Loading tables...")

    print("Loading {}".format(exc_table))
    exc_synapse = pd.read_csv(exc_table, index_col=0)
    exc_synapse.set_index('synapse_id', inplace=True)

    print("Loading {}".format(cell_table))
    cells_df = pd.read_csv(cell_table, index_col=0)
    cells_df.set_index('pt_root_id', inplace=True)

    print("Pre-loaded all tables")
    print("_"*80)

    # Loop through the directories in dir and generate metrics for each
    for subdir in os.listdir(_dir):
        if not os.path.isdir(_dir + '/' + subdir):
            continue

        print("Generating metrics for {}".format(_dir + '/' + subdir))
        get_metrics(_dir + '/' + subdir, exc_synapse, cells_df, _ignore_burst=_ignore_burst)
        print("Finished generating metrics for {}".format(_dir + '/' + subdir))
        print("_"*80)


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='Analyze ground truth metrics')
    parser.add_argument('--dir', type=str, help='Directory containing ground truth synapse table')
    parser.add_argument('--exc_synapse', type=str, default='synapses_w_ids.csv', help='Synapse table to use for excitatory synapses')
    parser.add_argument('--cell_table', type=str, default='cells_no_repeats.csv', help='Cell table to use')
    parser.add_argument('--ignore_burst', type=bool, default=False, help='Ignore burst cells')

    args = parser.parse_args()

    main(args.dir, args.exc_synapse, args.cell_table, args.ignore_burst)
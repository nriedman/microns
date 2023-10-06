### A script to carry out the entire pipeline to query additional synapses for a given sample size and calculate the precision and recall of the model.

import pandas as pd
import numpy as np
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import utils as ut
from caveclient import CAVEclient
import nglui
from tqdm import tqdm
import os
import time


"""
STEP 1: Query the synapses for the ground truth the model
STEP 2: Clean and align the data
STEP 3: Generate the mst's for the ground truth and the model and extract the paths
STEP 4: Calculate the precision and recall
"""

def sample_cell_ids(synapse_counts, sample_count=1000, synapse_min=6):
    cell_ids_sample = np.array(synapse_counts.sample(sample_count).index.values)
    count = 0
    count_limit = 1000

    while min(synapse_counts.loc[cell_ids_sample].values) < synapse_min and count < count_limit:
        print(min(synapse_counts.loc[cell_ids_sample]))
        cell_ids_sample = cell_ids_sample[synapse_counts.loc[cell_ids_sample] >= synapse_min]
        cell_ids_sample = np.append(cell_ids_sample, synapse_counts.sample(sample_count - len(cell_ids_sample)).index.values)
        count += 1

    if count == count_limit:
        print(f'Could not find a cell with at least {synapse_min} synapses')
    else:
        print(f'Found {sample_count} cells with at least {synapse_min} synapses:')
        print(min(synapse_counts.loc[cell_ids_sample].values))

    return cell_ids_sample


def query_ground_truth(client, synapse_df, synapse_count=200, sample_size=100):
    wanted_columns = ['pre_pt_root_id', 'post_pt_root_id', 'size', 'ctr_pt_position']
    ground_truth_df = pd.DataFrame(columns=wanted_columns + ['ctr_pt_position_x', 'ctr_pt_position_y', 'ctr_pt_position_z'])
    participant_ids = []
    under_sampled = []

    count = 0
    count_limit = 100

    while len(participant_ids) < sample_size and count < count_limit:
        new_id_pool = synapse_df[~synapse_df['post_pt_root_id'].isin(participant_ids) & ~synapse_df['post_pt_root_id'].isin(under_sampled)]
        new_cell_ids = sample_cell_ids(new_id_pool.groupby('post_pt_root_id').size(), sample_count=sample_size - len(participant_ids))
        print(f'Num new cell ids: {len(new_cell_ids)}')

        if synapse_count == 0:
            participant_ids = new_cell_ids
            break
        elif synapse_count < 0:
            participant_ids = new_cell_ids
            new_synapses = client.materialize.query_table('synapses_pni_2',
                                                            filter_in_dict={'post_pt_root_id': new_cell_ids},
                                                            select_columns=wanted_columns,
                                                            split_positions=True)
            ground_truth_df = pd.concat([ground_truth_df, new_synapses], ignore_index=True)
            break

        for i, cell_id in enumerate(new_cell_ids):
            new_synapses = client.materialize.query_table('synapses_pni_2',
                                                            filter_equal_dict={'post_pt_root_id': cell_id},
                                                            select_columns=wanted_columns,
                                                            split_positions=True,
                                                            limit=synapse_count)
            
            if new_synapses.shape[0] < synapse_count:
                print(f'Cell {cell_id} only has {new_synapses.shape[0]} synapses')
                under_sampled.append(cell_id)
            else:
                participant_ids.append(cell_id)
                ground_truth_df = pd.concat([ground_truth_df, new_synapses], ignore_index=True)
            print(f'Finished {i+1} of {len(new_cell_ids)}', end='\r')
        print(f'Num unique cells: {ground_truth_df["post_pt_root_id"].nunique()}')
        count += 1
    
    if count == count_limit:
        print(f'Could not find {sample_size} cells')
    else:
        print(f'Found {sample_size} cells')
    
    print(f'Finished with {len(participant_ids)} cells')
    return ground_truth_df.drop(columns='ctr_pt_position'), participant_ids


def clean_data(ground_truth_df, slimmed_synapse_df, participant_ids):
    # Rename position columns
    name_map = {'ctr_pt_position_x': 'ctr_pt_x', 'ctr_pt_position_y': 'ctr_pt_y', 'ctr_pt_position_z': 'ctr_pt_z'}
    ground_truth_df = ground_truth_df.rename(columns=name_map)

    # Convert position columns to microns
    ground_truth_df[['ctr_pt_x', 'ctr_pt_y']] *= 4 / 1000
    ground_truth_df['ctr_pt_z'] *= 40 / 1000

    # Merge the synapse dataframes
    print("Merging synapse dataframes")
    ground_truth_df['synapse_id'] = np.nan
    ground_truth_df['cell_type_pre'] = 'Unknown'
    for i in tqdm(range(ground_truth_df.shape[0])):
        row = ground_truth_df.iloc[i]
        matching_rows = slimmed_synapse_df[(slimmed_synapse_df['pre_pt_root_id'] == row['pre_pt_root_id']) & (slimmed_synapse_df['post_pt_root_id'] == row['post_pt_root_id']) & (slimmed_synapse_df['ctr_pt_x'] == row['ctr_pt_x']) & (slimmed_synapse_df['ctr_pt_y'] == row['ctr_pt_y']) & (slimmed_synapse_df['ctr_pt_z'] == row['ctr_pt_z']) & (slimmed_synapse_df['size'] == row['size'])]
        if matching_rows.shape[0] == 1:
            ground_truth_df.at[i, 'synapse_id'] = matching_rows['synapse_id'].values[0]
            ground_truth_df.at[i, 'cell_type_pre'] = matching_rows['cell_type_pre'].values[0]
        elif matching_rows.shape[0] > 1:
            print('Error: multiple matching rows')
            print(matching_rows)
        else:
            ground_truth_df.at[i, 'synapse_id'] = -(i + 2)
    
    # Add the remaining rows from slimmed_synapse_df
    print("Adding remaining excitatory rows")
    for cell_id in tqdm(participant_ids):
        matching_rows = slimmed_synapse_df[(slimmed_synapse_df['post_pt_root_id'] == cell_id) & (~slimmed_synapse_df['synapse_id'].isin(ground_truth_df['synapse_id']))]
        ground_truth_df = pd.concat([ground_truth_df, matching_rows[ground_truth_df.columns]], ignore_index=True)
    
    return ground_truth_df


def prepare_ground_truth(client, synapse_df, dir, sample_size, synapse_count, version):
    # Create output directory
    output = dir + '/' + version + '_p' + str(sample_size) + '_s' + str(synapse_count) + '_a0'
    attempt = 1
    const_len = len(output) - 1
    while os.path.exists(output):
        output = output[:const_len] + str(attempt)
        attempt += 1
    os.mkdir(output)

    ground_truth_df, participant_ids = query_ground_truth(client, synapse_df, synapse_count=synapse_count, sample_size=sample_size)

    # Check that the number of cells is correct
    assert len(participant_ids) == sample_size

    # Save the participant ids
    pd.to_pickle(participant_ids, output + '/participant_ids.pkl')

    # Clean and align the data
    slimmed_exc_df = synapse_df[synapse_df['post_pt_root_id'].isin(participant_ids)]
    ground_truth_df = clean_data(ground_truth_df, slimmed_exc_df, participant_ids)
    print("Finished cleaning data")

    # Save the ground truth
    print("Saving ground truth synapses to {}".format(output + '/ground_truth_synapses.csv'))
    ground_truth_df.to_csv(output + '/ground_truth_synapses.csv')

    print("Finished preparing ground truth")
    print("-"*80)


def main(dir, sample_size=100, synapse_counts=[200], repititions=1, version='v343', exc_synapse='synapses_w_ids.csv'):
    # Initialize client
    client = CAVEclient('minnie65_public_' + version)
    queries = 0

    # Query the synapses for the ground truth
    synapse_df = pd.read_csv(exc_synapse, index_col=0)

    print(synapse_counts)
    for synapse_count in synapse_counts:
        for trial in range(repititions):
            if queries > 2000:
                print("Waiting for 10 minutes...")
                time.sleep(600)
                queries = 0

            print(f'Preparing ground truth for {synapse_count} synapses, trial {trial+1}/{repititions}')
            prepare_ground_truth(client, synapse_df, dir, sample_size, synapse_count, version)
            queries += sample_size
    
    print("Finished all trials")
    

if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='Evaluate the model')
    parser.add_argument('--dir', type=str, default=os.getcwd(), help='Directory to save the output to')
    parser.add_argument('--sample_size', type=int, default=100, help='Number of cells to sample')
    parser.add_argument('--synapse_count', type=int, default=200, nargs='+', help='Number of synapses to query per cell. Default is 200')
    parser.add_argument('--repititions', type=int, default=1, help='Number of times to repeat the experiment for each synapse_count. Default is 1')
    parser.add_argument('--version', type=str, default='v343', help='Data version to use')
    parser.add_argument('--exc_synapse', type=str, default='synapses_w_ids.csv', help='Synapse table to use for excitatory synapses')

    args = parser.parse_args()

    main(args.dir, args.sample_size, args.synapse_count, args.repititions, args.version, args.exc_synapse)
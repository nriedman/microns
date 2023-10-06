import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import utils as ut
from tqdm import tqdm
import pylcs as LCS
import concurrent.futures
import argparse
import os


def shuffle_sequences(strings):
    return np.array([''.join(np.random.permutation(list(string))) for string in strings])


# Function to calculate LCS sizes for a range of indices (GPT round 2)
def calculate_lcs_range(start_idx, end_idx, sorted_strings):
    lcs_dict = {}
    for i in tqdm(range(start_idx, end_idx), desc="Processing"):
        str1 = sorted_strings[i]
        length1 = len(str1)
        lcs_dict.setdefault(length1, [])

        lcs_comps = np.array(LCS.lcs_string_of_list(str1, sorted_strings[i+1:]))

        max_lcs = max(lcs_comps)
        lcs_seq_pairs = sorted_strings[i+1:][lcs_comps == max_lcs]

        # Add the sequence, max LCS, and similar sequences to the dictionary
        lcs_dict[length1].append((str1, max_lcs, lcs_seq_pairs))
    return lcs_dict


def calculate_lcs(sorted_strings, num_processes = 4):
    process_start = 0
    process_end = len(sorted_strings)-1
    process_range = process_end - process_start
    
    # Split up range into chunks for each process
    chunk_size = process_range // num_processes
    ranges = [(process_start + i * chunk_size, process_start + (i + 1) * chunk_size) for i in range(num_processes - 1)] + [(process_start + (num_processes - 1) * chunk_size, process_end)]
    
    # Parralelize LCS calculation
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [executor.submit(calculate_lcs_range, start, end, sorted_strings) for start, end in ranges]

    # Wait for all futures to complete
    concurrent.futures.wait(futures)

    # Combine results
    lcs_dict = {}
    for future in futures:
        for length, lcs_list in future.result().items():
            lcs_dict.setdefault(length, [])
            lcs_dict[length].extend(lcs_list)
    
    return lcs_dict



def main(input_file, output_dir, length, num_trials, num_processes=4):
    # Load the data
    # A list of encodded sequence strings
    print("Loading data...")
    all_sequences = pd.read_pickle(input_file)
    all_sequences.sort(key=len)
    print("Loaded {} sequences.".format(len(all_sequences)))
    
    sequences = [seq for seq in all_sequences if len(seq) >= length]
    print("Filtered to {} sequences of length {} or greater.".format(len(sequences), length))
    print("-"*50)

    print('Beginning trials...')
    for i in range(num_trials):
        print("Trial {} of {}...".format(i+1, num_trials))
        print("-"*50)

        # Shuffle the sequences
        print("Shuffling sequences...")
        shuffled_sequences = shuffle_sequences(sequences)

        # Calculate LCS values
        print("Calculating LCS values...")
        new_data = calculate_lcs(shuffled_sequences, num_processes=num_processes)

        # Save the data
        print("Saving data...")

        output = os.path.join(output_dir, 'l_{}_t_1_a_0'.format(length))
        while os.path.exists(output + '.pkl'):
            output = output[:-1] + str(int(output[-1]) + 1)
        output += '.pkl'

        pd.to_pickle(new_data, output)
        print("Saved data to {}".format(output))
        print("-"*50)
    
    print("Done.")


if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser(description='Generate the null model data.')
    parser.add_argument('--input', type=str, help='The path to the encoded sequences file (pickle).')
    parser.add_argument('--output', type=str, help='The path to the output directory.')
    parser.add_argument('--length', type=int, help='The minimum length of the sequences to consider.')
    parser.add_argument('--num_trials', type=int, help='The number of trials to run.')
    parser.add_argument('--num_processes', type=int, default=4, help='The number of processes to use.')

    args = parser.parse_args()

    main(args.input, args.output, args.length, args.num_trials, args.num_processes)
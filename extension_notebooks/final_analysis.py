# Path: final_analysis.py

from tqdm import tqdm
import numpy as np
import pandas as pd
import pylcs as LCS
import argparse
import os
import concurrent.futures


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

def main(input_file, output_dir, _process_start, _process_end, num_processes):
    # Read in the encoded sequences (pickle file)
    encoded_sequences = pd.read_pickle(input_file)

    # Sort the sequences by length
    sorted_strings = sorted(encoded_sequences, key=len)

    # Designate start and end indexes for the range of indices to process
    process_start = int(_process_start * (len(sorted_strings)-1))
    process_end = int(_process_end * (len(sorted_strings)-1))
    process_range = process_end - process_start

    print('Processing sequences {} to {} ({} total)'.format(process_start, process_end, process_range))

    # Split the work into chunks for parallel processing
    print("Splitting work into {} chunks...".format(num_processes))
    chunk_size = process_range // num_processes
    ranges = [(process_start + i * chunk_size, process_start + (i + 1) * chunk_size) for i in range(num_processes - 1)] + [(process_start + (num_processes - 1) * chunk_size, process_end)]

    print("Split work into the following ranges:")
    print(ranges)

    # Parallelize the LCS calculations using ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [executor.submit(calculate_lcs_range, start, end, np.array(sorted_strings)) for start, end in ranges]

    # Wait for all futures to complete
    concurrent.futures.wait(futures)

    print("Analysis complete. Writing to file...")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save the results from each process to a separate pickle file
    for i, future in enumerate(futures):
        pd.to_pickle(future.result(), output_dir + "/s_{}_e_{}_c_{}_.pkl".format(int(_process_start * 100), int(_process_end*100), i))
    
    print("Done.")


if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate LCS sizes for a range of indices')
    parser.add_argument('--input', type=str, help='Path to the encoded sequence list file (pickle)')
    parser.add_argument('--output', type=str, help='Path to the output directory')
    parser.add_argument('--start', type=float, help='Start index for the range of indices to process (float: percent of all sequences)', default=0.0)
    parser.add_argument('--end', type=float, help='End index for the range of indices to process (float: percent of all sequences)', default=1.0)
    parser.add_argument('--processes', type=int, help='Number of processes to use', default=4)
    
    args = parser.parse_args()
    main(args.input, args.output, args.start, args.end, args.processes)
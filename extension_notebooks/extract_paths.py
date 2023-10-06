
import pandas as pd
import numpy as np
import utils as ut
import os
import argparse
import math
import string


def map_to_printable_string(number, length=3):
    if 0 <= number <= 56208:
        printable_chars = string.printable[:-6]  # Exclude non-printable characters
        base = len(printable_chars)
        
        encoded = ""
        while number > 0:
            number, index = divmod(number, base)
            encoded = printable_chars[index] + encoded
        
        # Pad with leading zeros if needed
        padding = length - len(encoded)
        encoded = printable_chars[0] * padding + encoded
        
        return encoded
    else:
        raise ValueError("Number must be in the range 0 to 56208")


def get_printable(encoded_str, encoding_map):
    return "".join([encoding_map[c] for c in encoded_str])


def encode_sequences(sequences_raw, cells_df):
    # Printable unicode blocks:
    # 19968 - 40959
    coding_block_1 = np.arange(19968, 40960, 1)

    # 44032 - 55215
    coding_block_2 = np.arange(44032, 55216, 1)

    # 131072 - 173791
    coding_block_3 = np.arange(131072, 173792, 1)
    code_points = np.concatenate((coding_block_1, coding_block_2, coding_block_3))

    cell_ids = cells_df.index.values

    # Create a dictionary mapping cell ids to code points
    pt_root_id_to_char = {pt_root_id: chr(code_points[i]) for i, pt_root_id in enumerate(cell_ids)}
    char_to_pt_root_id = {v: k for k, v in pt_root_id_to_char.items()}

    # Create a dictinoary mapping chars to printable codons
    char_to_codon = {char: map_to_printable_string(i) for i, char in enumerate(pt_root_id_to_char.values())}

    # Convert raw sequences to encoded sequences
    sequences_encoded = {}
    for cell_id, sequences in sequences_raw.items():
        sequences_encoded[cell_id] = []
        for sequence in sequences:
            encoded = "".join([pt_root_id_to_char[pt_root_id] for pt_root_id in sequence])
            sequences_encoded[cell_id].append(encoded)
    
    return sequences_encoded, pt_root_id_to_char, char_to_pt_root_id, char_to_codon


def main(dataset, cell_table, dir):
    # # Load dataset
    # print("Loading datasets...")

    # print("Loading {}".format(dataset))
    # syn_df = pd.read_csv(dataset, index_col=0)
    # syn_df.set_index('synapse_id', inplace=True)

    print("Loading {}".format(cell_table))
    cell_df = pd.read_csv(cell_table, index_col=0)
    cell_df.set_index('pt_root_id', inplace=True)

    # print("Loaded all tables")
    # print("_"*80)

    # # Generate msts
    # print("Generating MSTs...")
    # msts, _ = ut.generate_msts(syn_df, cell_df)

    # # Extract the paths from the msts
    # print("Extracting unencoded sequences...")
    # sequences_raw = {}
    # for i, mst in enumerate(msts):
    #     cell_id = mst.graph['cell_id']
    #     sequences_raw[cell_id] = []

    #     print(cell_id, f'{i+1}/{len(msts)}', end='\r')
    #     new_sequences = ut.get_paths(mst, root=-1, ignore_burst=True, return_seq=True)[1]
    #     for seq in new_sequences:
    #         if len(seq) > 2:
    #             sequences_raw[cell_id].append(seq)
    
    # print("Sequences extracted. Saving raw sequences...")
    # # Save the raw sequences
    # if not os.path.exists(dir):
    #     os.makedirs(dir)
    
    # pd.to_pickle(sequences_raw, dir + '/sequences_raw.pkl')

    sequences_raw = pd.read_pickle(dir + '/sequences_raw.pkl')

    print("Raw sequences saved. Encoding sequences as strings...")

    # Encode the sequences
    encoded_sequences, pt_root_id_to_char, char_to_pt_root_id, char_to_codon = encode_sequences(sequences_raw, cell_df)

    print("Sequences encoded. Saving encoded sequences...")
    # Save the encoded sequences as a dictionary
    pd.to_pickle(encoded_sequences, dir + '/sequences_encoded_dict.pkl')
    print("Saved encoded sequence dict to {}".format(dir + '/sequences_encoded_dict.pkl'))

    # Save the maps
    pd.to_pickle(pt_root_id_to_char, dir + '/pt_root_id_to_char.pkl')
    print("Saved pt_root_id_to_char to {}".format(dir + '/pt_root_id_to_char.pkl'))

    pd.to_pickle(char_to_pt_root_id, dir + '/char_to_pt_root_id.pkl')
    print("Saved char_to_pt_root_id to {}".format(dir + '/char_to_pt_root_id.pkl'))

    pd.to_pickle(char_to_codon, dir + '/char_to_codon.pkl')
    print("Saved char_to_codon to {}".format(dir + '/char_to_codon.pkl'))

    # Save a list of the encoded sequences
    encoded_sequences_list = []
    for cell_id, sequences in encoded_sequences.items():
        for sequence in sequences:
            encoded_sequences_list.append(sequence)
    
    pd.to_pickle(encoded_sequences_list, dir + '/sequences_encoded_list.pkl')
    print("Saved encoded sequence list to {}".format(dir + '/sequences_encoded_list.pkl'))


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='Extract paths from a given dataset')
    parser.add_argument('--syn_table', type=str, default='synapses_w_ids.csv', help='Synapse dataset to extract paths from')
    parser.add_argument('--cell_table', type=str, default='cells_no_repeats.csv', help='Dataset with cell info')
    parser.add_argument('--dir', type=str, default='sequences.csv', help='Output directory')
    parser.add_argument('--ignore_bursts', type=bool, default=False, help='Ignore burst sequences')

    args = parser.parse_args()

    main(args.syn_table, args.cell_table, args.dir)
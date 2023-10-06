import pandas as pd
import numpy as np
import argparse
import os


def merge_dicts(dict_list):
    result_dict = {}
    
    for dictionary in dict_list:
        for key, value in dictionary.items():
            if key in result_dict:
                result_dict[key].extend(value)
            else:
                result_dict[key] = value
    
    return result_dict


def save_pickle(data, markers, out_dir):
    filename = ''.join([str(m[0]) + '_' + str(m[1]) + '_' for m in markers])
    output = os.path.join(out_dir, filename + 'a_0')

    while os.path.exists(output + '.pkl'):
        output = output[:-1] + str(int(output[-1]) + 1)
    
    output += '.pkl'
    pd.to_pickle(data, output)
    print("Saved data to {}".format(output))


def main(dir, out_dir, spent, null=False):
    # Get the data files
    data_ls = []
    for file in os.listdir(dir):
        if file.endswith(".pkl"):
            data_ls.append(os.path.join(dir, file))

    data = []
    for data_file in data_ls:
        data.append(pd.read_pickle(data_file))
    print("Loaded files: {} ({} total)".format([os.path.basename(f) for f in data_ls], len(data_ls)))

    # Merge the dictionaries
    merged_data = merge_dicts(data)
    print("Merged {} files into a single dictionary.".format(len(data_ls)))

    if not null:
        # Get the start and end percentages from the file names, marked by underscores with keys 's' and 'e'
        starts, ends = [], []
        for path in data_ls:
            filename = os.path.basename(path)
            components = np.array(filename.split('_'))

            starts.append(int(components[np.where(components == 's')[0][0] + 1]))
            ends.append(int(components[np.where(components == 'e')[0][0] + 1]))

        # Save the data
        save_pickle(merged_data, zip(['s', 'e'], [min(starts), max(ends)]), out_dir)

    else:
        # Get the number of trials from the file names and add them together, marked by underscores with key 't'
        trial_sum = 0
        for path in data_ls:
            filename = os.path.basename(path)
            components = np.array(filename.split('_'))
            
            trial_sum += int(components[np.where(components == 't')[0][0] + 1])
        
        length = int(components[np.where(components == 'l')[0][0] + 1])

        # Save the data
        save_pickle(merged_data, zip(['l', 't'], [length, trial_sum]), out_dir)


    # Move the data files to the combined 'spent' directory
    for file in data_ls:
        os.rename(file, os.path.join(spent, os.path.basename(file)))
        print("Moved {} to {}".format(file, os.path.join(spent, os.path.basename(file))))


if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser(description='Merge the data from the given files into a single file.')
    parser.add_argument('--in_dir', type=str, help='The directory containing the data files to merge (pickle).')
    parser.add_argument('--out_dir', type=str, help='The path to the output directory.')
    parser.add_argument('--spent', type=str, help='The path to the spent directory, where the data files will be moved to after merging.')
    parser.add_argument('--null', type=bool, default=False, help='Whether or not the data is from the null model. If not, the data will be merged into a single file with the start and end percentages in the file name.')

    args = parser.parse_args()

    main(args.in_dir, args.out_dir, args.spent, args.null)
    
# Microns

Summer 2023 project analyzing the synaptic organization of the dendrite in the MICrONS Cortial MM^3 dataset.

sample_dataset:
- .csv files containing all the input synapses of three 23P cells
- The number in the filename corresponds to the excitatory input count percentile of that cell

cells_no_repeats.csv:
- Contains cell information on all the excitatory cells is in the entire MM^3 dataset
- Needed to generate msts

gen_msts.py:
- Python script including the function and dependencies for generating minimum spanning trees from a given synapse table

**For Generating Msts (gen_msts.py)**

The function that generates minimum spanning trees is called **generate_msts()**.

The arguments (described in detail in the script) are:
- synapses_df: A pandas dataframe of synapses
- cells_df: A pandas dataframe of cell info
- k: An integer, the number of nearest neighbors for each synapse to consider when generating the minimum spanning tree
- soma_k: An integer, the number of nearest neighbors for the soma to consider when generating the minimum spanning tree

It returns **two** objects:
- First: A list of minimum spanning trees, whose properties are described in detail in the script.
- Second: A list of the cell_ids that did not have enough synapses to make a minimum spanning tree (likely empty, often ignored)

**Notes on usage:**

- The synapse table that you pass into the function can include data for more than one post-cell. The function will return a list of msts, one for each post-cell in the synapse table.

- The synapse dataframe must be indexed by synapse_id, and the cell dataframe must be indexed by pt_root_id.

- The function needs each dataframe to include certain columns. An exact list is in the function description in the list.

Here is an example of how I would use the function:

https://github.com/nriedman/microns/blob/7d1233d3ac5eba714d28369980fcdc283de6c08c/mst_example.py

**For Plotting Msts (plot_msts.py)**

The plot_msts.py script has functions for visualizing msts in 3D, highlighting sequences.

The function to do so is called plot_mst_3d(). It takes 1 positional argument (G, an mst graph) and returns the paths for that graph. Along the way, it plots the graph in 3D, highlighting paths of synapses in color. I included it in case you want to visualize any of the minimum spanning trees you create.

plot_msts.py also includes the functions to get the paths from a given minimum spanning tree. I added functionality to ignore burst sequences, as well as an option to directly return sequences of cell_ids instead of paths of synapse ids. If you get sequences from a tree and would like to conver them into an easily readible format, I also uploaded the dictionaries I use that map pre-cell ids to single unicode characters and back again, as well as from the single unicode characters to a unique 3-character codon that is guaranteed to be legible.

- pre-cell id to char: pt_root_id_to_char.pkl

- char to pre-cell id: char_to_pt_root_id.pkl

- char to codon: char_to_codon.pkl

Let me know if you have any questions!


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

The function that generates minimum spanning trees is called **generate_msts()**.

The arguments (described in detail in the script) are:
- synapses_df: A pandas dataframe of synapses
- cells_df: A pandas dataframe of cell info
- k: An integer, the number of nearest neighbors for each synapse to consider when generating the minimum spanning tree
- soma_k: An integer, the number of nearest neighbors for the soma to consider when generating the minimum spanning tree

It returns a list of minimum spanning trees, whose properties are described in detail in the script.

**Notes on usage:**

- The synapse table that you pass into the function can include data for more than one post-cell. The function will return a list of msts, one for each post-cell in the synapse table.

- The synapse dataframe must be indexed by synapse_id, and the cell dataframe must be indexed by pt_root_id.

- The function needs each dataframe to include certain columns. An exact list is in the function description in the list.

Here is an example of how I would use the function:




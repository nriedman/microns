import gen_msts as gm

# Load in desired synapse table, setting index_col to 0 because I always forget to ignore the index when I save csvs
synapse_df = pd.read_csv('23P_corner_95.csv', index_col=0)

# Set the index to synapse_id
synapse_df.set_index('synapse_id', inplace=True)

# Repeat for the cell_table
cells_df = pd.read_csv('cells_no_repeats.csv', index_col=0)
cells_df.set_index('pt_root_id', inplace=True)

# Optional: Slim the cells_df to only include information on the cells present in synapse_df
cells_df = cells_df[cells_df['pt_root_id'].isin(synapse_df['post_pt_root_id'].unique())]

# Generate the msts for the given data (I chose not to provide k, soma_k and instead allowed them to be default)
msts, _ = gm.generate_msts(synapse_df, cells_df)

# I now have a list of minimum spanning trees (msts)

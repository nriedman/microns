import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.neighbors import NearestNeighbors
from copy import deepcopy
import pylcs as LCS
from tqdm import tqdm


"""
Sequence algorithm:
All about finding sequences from the minimum spanning tree of the synapse space (AKA dendrite)
"""

def traverse_branch(G, start_node, root_node, ignore_burst=False):
    path = []
    node = start_node
    safeguard = root_node

    while G.degree(node) == 2:
        # If we know the cell type, it's excitatory, so add it to the path
        node_ct = G.nodes[node]['cell_type']
        cur_cell_id = G.nodes[node]['pre_cell_id']
        prev_cell_id = G.nodes[safeguard]['pre_cell_id']

        if node_ct != 'Unknown':
            if ignore_burst:
                if safeguard == root_node or cur_cell_id != prev_cell_id:
                    path.append(node)
            else:
                path.append(node)
        
        # Make sure we don't go backwards along the path
        neighbors = np.array(list(G.neighbors(node)))
        next_node = neighbors[neighbors != safeguard][0]

        # March forward a step
        safeguard = node
        node = next_node
    
    # Add the last node to the path
    node_ct = G.nodes[node]['cell_type']
    cur_cell_id = G.nodes[node]['pre_cell_id']
    prev_cell_id = G.nodes[safeguard]['pre_cell_id']

    if node_ct != 'Unknown':
        if ignore_burst:
            if safeguard == root_node or cur_cell_id != prev_cell_id:
                path.append(node)
        else:
            path.append(node)

    return path, node, safeguard

# Given a graph and a root node, return a list of lists of nodes, where each list is a sequence of nodes such that each node is connected to the next node in the list and the last node in the list has degree > 2
def get_paths(G, root, entrance_node=None, ignore_burst=False, return_seq=False):
    paths = []
    sequences = []

    root_neighbors = np.array(list(G.neighbors(root)))
    new_branches = root_neighbors[~np.isin(root_neighbors, entrance_node)]

    for node in new_branches:
        path, end_node, visited_node = traverse_branch(G, node, root, ignore_burst=ignore_burst)
        paths.append(path)
        sequences.append([G.nodes[n]['pre_cell_id'] for n in path])

        if G.degree(end_node) > 2:
            if return_seq:
                new_paths, new_sequences = get_paths(G, end_node, visited_node, ignore_burst=ignore_burst, return_seq=return_seq)
                sequences.extend(new_sequences)
            else:
                new_paths = get_paths(G, end_node, visited_node, ignore_burst=ignore_burst, return_seq=return_seq)
            
            paths.extend(new_paths)
    
    if return_seq:
        return paths, sequences
    else:
        return paths


def get_color_map(paths, root):
    colors = ['#ff0000', '#00ff00', '#0000ff', '#ffff00', '#00ffff', '#ff00ff', '#ff8000', '#8000ff', '#00ff80', '#ff0080']
    color_map = {node: colors[i % len(colors)] for i, path in enumerate(paths) for node in path}
    color_map[root] = '#a0a0a0'
    return color_map


# Plot the synapes of the minimum spanning tree in 3D and connect them with lines according to the minimum spanning tree
def plot_mst_3d(G, root=-1, paths=None, ignore_burst=False):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    if paths is None:
        paths = get_paths(G, root, ignore_burst=ignore_burst)
    
    color_map = get_color_map(paths, root)

    for node in G.nodes():
        if node in color_map:
            color = color_map[node]
        else:
            color = '#999999'
        coords = G.nodes[node]['pos']
        ax.scatter(coords[0], coords[1], coords[2], c=color, s=10)
    
    # Connect the synapses with lines according to the edges of the minimum spanning tree
    for edge in G.edges():
        node1, node2 = G.nodes[edge[0]]['pos'], G.nodes[edge[1]]['pos']
        ax.plot([node1[0], node2[0]],
                [node1[1], node2[1]],
                [node1[2], node2[2]],
                c='k')
    
    # Plot the soma
    soma_pos = G.nodes[root]['pos']
    ax.scatter(soma_pos[0], soma_pos[1], soma_pos[2], c='r', s=100)

    ax.set_xlabel('z')
    ax.set_ylabel('y')
    ax.set_zlabel('x')

    plt.tight_layout()
    plt.show()

    return paths


"""
Generating the minimum spanning tree:
"""

# Given a graph G, find all its connected components, and for each component not connected to the root, add an edge between the soma and the closest leaf node in the component to the soma
def connect_disjoint_branches(G, soma_node=-1):
    G = deepcopy(G)
    G.remove_edges_from(list(nx.selfloop_edges(G)))

    components = list(nx.connected_components(G))
    root_component = [c for c in components if soma_node in c][0]
    components.remove(root_component)

    soma_pos = G.nodes[soma_node]['pos']

    for component in components:
        component = list(component)
        leaf_nodes = [node for node in component if G.degree(node) == 1]
        distances = [np.linalg.norm(G.nodes[leaf_node]['pos'] - soma_pos) for leaf_node in leaf_nodes]

        closest_leaf_node = leaf_nodes[np.argmin(distances)]
        G.add_edge(soma_node, closest_leaf_node)
    
    return G



def generate_msts(synapses_df, cells_df, k=6, soma_k=10):
    """Generate a minimum spanning tree for each cell in the dataset.
    
    Parameters
    ----------
    synapses_df : pandas.core.frame.DataFrame
        A dataframe of synapse data.
        Must include columns for ctr_pt_x, ctr_pt_y, ctr_pt_z, post_pt_root_id, pre_pt_root_id, and pre_cell_type.
        Index should be the synapse id.
    cells_df : pandas.core.frame.DataFrame
        A dataframe of cell data, including columns for pt_root_id, pt_x, pt_y, pt_z, and cell_type.
        Index should be the pt_root_id.
    k : int
        The number of nearest neighbors to use when constructing the minimum spanning tree.
    soma_k : int
        The number of nearest neighbors to use when constructing the minimum spanning tree for the soma.
    Returns
    -------
    msts : list
        A list of minimum spanning trees.
        Each tree is a networkx.Graph object. The nodes are synapse ids, and the edge weights are the distances between synapses.
        Graph attributes:
            G.graph['cell_id'] : int
                The cell id of the post-synaptic cell the synapses belong to.
            G.graph['cell_type'] : str
                The cell type of the post-synaptic cell the synapses belong to.
        Node attributes:
            G.nodes[syn_id]['pos'] : numpy.ndarray
                The xyz coordinates of the synapse.
            G.nodes[syn_id]['cell_type'] : str
                The cell type of the pre-synaptic cell the synapse belongs to.
            G.nodes[syn_id]['pre_cell_id'] : int
                The cell id of the pre-synaptic cell the synapse belongs to.
        Soma node attributes:
            G.nodes[-1]['pos'] : numpy.ndarray
                The xyz coordinates of the soma.
    """
    synapses_grouped = synapses_df.groupby('post_pt_root_id')
    msts = []
    too_sparse = []
    count = 0
    total = len(synapses_grouped)
    for cell_id, syn_group in synapses_grouped:
        # If there are less than k synapses, skip this cell
        group_size = syn_group.shape[0]
        if group_size < k:
            too_sparse.append(cell_id)
            continue

        # Keep relevant rows of synapse table
        synapses = syn_group[['ctr_pt_x', 'ctr_pt_y', 'ctr_pt_z']]

        # Get the soma location for the cell
        cell_info = cells_df.loc[cell_id]
        soma_xyz = np.array(cell_info[['pt_x', 'pt_y', 'pt_z']].values)
        soma_xyz = np.matmul(soma_xyz, np.diag([4/1000, 4/1000, 40/1000]))

        # Add the soma location to the synapse table
        soma_df = pd.DataFrame(soma_xyz).T
        soma_df.columns = ['ctr_pt_x', 'ctr_pt_y', 'ctr_pt_z']
        soma_df.index = [-1]
        synapses_w_soma = pd.concat([synapses, soma_df])

        # Create a kdtree from the synapse locations
        kd_tree = NearestNeighbors(n_neighbors=k, algorithm='kd_tree').fit(synapses_w_soma.values)

        # Get the k nearest neighbors for each synapse and the soma
        distances, indices = kd_tree.kneighbors(synapses.values)
        if group_size >= soma_k:
            soma_distances, soma_indices = kd_tree.kneighbors(soma_xyz.reshape(1, -1), n_neighbors=soma_k)
        else:
            soma_distances, soma_indices = kd_tree.kneighbors(soma_xyz.reshape(1, -1), n_neighbors=group_size)

        # Subtract the "radius" of the soma from the soma distances
        soma_radius = soma_distances[0][1]
        soma_distances = soma_distances - soma_radius 

        # Create a graph from the synapse group
        nodes = list(synapses.index.values)
        G = nx.Graph(cell_id=cell_id, cell_type=cells_df.loc[cell_id, 'cell_type'])
        for node in nodes:
            node_ct = syn_group.loc[node, 'cell_type_pre']
            node_pre_id = syn_group.loc[node, 'pre_pt_root_id']
            G.add_node(node, pos=synapses.loc[node, ['ctr_pt_x', 'ctr_pt_y', 'ctr_pt_z']].values,
                             cell_type=node_ct,
                             pre_cell_id=node_pre_id)
        nodes.append(-1)
        G.add_node(-1, pos=soma_xyz, pre_cell_id=-1)

        # Add edges according to the kdtree
        for i in range(len(indices)):
            syn_id = nodes[i]
            for j in range(len(indices[i])):
                if i != indices[i][j]:
                    G.add_edge(syn_id, nodes[indices[i][j]], weight=distances[i][j])
        
        # Add edges from the soma to its nearest neighbors, corrected for the radius of the soma
        # Make sure not to add an edge from the soma to itself
        for l in range(1,len(soma_indices[0])):
            G.add_edge(-1, nodes[soma_indices[0][l]], weight=soma_distances[0][l])

        # Get the minimum spanning tree
        mst = nx.minimum_spanning_tree(G)

        # Make the graph fully connected
        mst = connect_disjoint_branches(mst)

        # Add the graph to the list of msts
        msts.append(mst)
        
        count += 1
        print(f'{cell_id}, {count}/{total}, {count/total}: {mst.number_of_nodes()} nodes, {mst.number_of_edges()} edges', end='\r')
    
    print('Done!')
    return msts, too_sparse

"""
For parallelization of extracting sequence overlap data.
"""
def calc_lcs(short_sequence, long_sequences):
    # max_lcs = 0
    # for long_sequence in long_sequences:
    #     lcs = LCS.lcs_string_length(short_sequence, long_sequence)
    #     if lcs == len(short_sequence):
    #         return lcs
    #     max_lcs = max(max_lcs, lcs)
    return max(LCS.lcs_string_of_list(short_sequence, long_sequences))

def calculate_distribution(sequence_lengths, sequences_sorted, k):
    kth_distribution = []
    cur_sequences = np.array(sequences_sorted[sequence_lengths[k]:sequence_lengths[k+1]] if k < len(sequence_lengths) - 1 else sequences_sorted[sequence_lengths[k]:])
    upper_sequences = np.array(sequences_sorted[sequence_lengths[k]:])

    for i in tqdm(range(len(cur_sequences))):
        # print(k, str(i * 100 / len(cur_sequences))[:4] + '%')
        kth_distribution.append(calc_lcs(cur_sequences[i], np.delete(upper_sequences, i)))
    return k, kth_distribution


# Function to calculate LCS sizes for a range of indices (GPT round 2)
def calculate_lcs_range(start_idx, end_idx, sorted_strings, lcs_sizes_by_length):
    for i in tqdm(range(start_idx, end_idx), desc="Processing"):
        str1 = sorted_strings[i]
        length1 = len(str1)
        lcs_sizes = lcs_sizes_by_length.setdefault(length1, [])
        lcs_sizes.append(LCS.lcs_string_of_list(str1, sorted_strings[i+1:]))

        # for j in range(i + 1, len(sorted_strings)):   # Vectorize this
        #     str2 = sorted_strings[j]
        #     lcs_size = LCS.lcs_string_length(str1, str2) # calc_lcs(str1, str2)
        #     lcs_sizes.append(lcs_size)


def get_printable(encoded_str, encoding_map):
    return "".join(['<' + encoding_map[c] + '>' for c in encoded_str])

if __name__ == '__main__':
    pass
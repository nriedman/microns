from copy import deepcopy
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.neighbors import NearestNeighbors


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
            G.nodes[-1]['pre_cell_id'] : int
                The pre-cell id of the soma (always -1, since soma has no pre-cell id).
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
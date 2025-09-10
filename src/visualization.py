import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# Given a graph and a root node, return a list of nodes such that each node is connected to the next node in the list and the last node in the list has degree > 2
def traverse_branch(G, start_node, root_node, ignore_burst=False, return_seq=False):
    path = []
    node = start_node
    safeguard = root_node

    while G.degree(node) == 2:
        # If we know the cell type, it's excitatory, so add it to the path
        node_ct = G.nodes[node]['cell_type']
        cur_cell_id = G.nodes[node]['pre_cell_id']
        prev_cell_id = G.nodes[safeguard]['pre_cell_id']

        if node_ct != 'Unknown':
            # If we're ignoring bursts, only add the node if it's a new cell
            if ignore_burst:
                if safeguard == root_node or cur_cell_id != prev_cell_id:
                    # If we're returning a sequence of cell ids, append the cell id, otherwise append the node
                    if return_seq:
                        path.append(cur_cell_id)
                    else:
                        path.append(node)
            else:
                if return_seq:
                    path.append(cur_cell_id)
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

    # If we know the cell type, it's excitatory, so add it to the path
    if node_ct != 'Unknown':
        # If we're ignoring bursts, only add the node if it's a new cell
        if ignore_burst:
            if safeguard == root_node or cur_cell_id != prev_cell_id:
                # If we're returning a sequence of cell ids, append the cell id, otherwise append the node
                if return_seq:
                    path.append(cur_cell_id)
                else:
                    path.append(node)
        else:
            if return_seq:
                path.append(cur_cell_id)
            else:
                path.append(node)

    return path, node, safeguard

# Given a graph and a root node, return a list of lists of nodes, where each list is a sequence of nodes such that each node is connected to the next node in the list and the last node in the list has degree > 2
def get_paths(G, root, entrance_node=None, ignore_burst=False, return_seq=False):
    paths = []

    # Identify the children of the root node, excluding the entrance node (prevent backtracking)
    root_neighbors = np.array(list(G.neighbors(root)))
    new_branches = root_neighbors[~np.isin(root_neighbors, entrance_node)]

    # Traverse each branch and add the paths to the list
    for node in new_branches:
        path, end_node, visited_node = traverse_branch(G, node, root, ignore_burst=ignore_burst, return_seq=return_seq)
        paths.append(path)

        if G.degree(end_node) > 2:
            # If the end node has degree > 2, we have a new branch, so recurse
            new_paths = get_paths(G, end_node, visited_node, ignore_burst=ignore_burst, return_seq=return_seq)
            paths.extend(new_paths)
    
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

    # If no paths are provided, get them from the graph
    if paths is None:
        paths = get_paths(G, root, ignore_burst=ignore_burst)
    
    color_map = get_color_map(paths, root)

    # Plot each synapse as a point
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
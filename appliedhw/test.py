import random

def BA(N, m):
    """
        Inputs:
            N = size of the network
            m = number of edges emanating from each newly added node

        Returns:
            Network = A dictionary with the adjacency list of the network
    """

    # How to construct the adjecency list

    # Initialize adjacency list (library)
    Network = {v: [] for v in range(N)} # creates list for each node v until N-1

    # Initial nucleus: fully connected clique of size m+1
    n0 = m + 1
    for v1 in range(n0):
        for v2 in range(v1 + 1 , n0):  # Avoid self-loops and duplicate edges
            Network[v1].append(v2)
            Network[v2].append(v1)

    # Auxiliary list for preferential attachment
    # Start by populating it with the initial nucleus nodes
    aux_list = []
    for v in range(n0):
        aux_list.extend([v] * len(Network[v]))

    # Add new nodes to the network
    for new_node in range(n0, N):
        connected_nodes = []
        while len(connected_nodes) < m:
            # Select a node based on preferential attachment
            selected_node = random.choice(aux_list)
            if selected_node not in connected_nodes:
                # Connect new_node to selected_node
                connected_nodes.append(selected_node)
                Network[new_node].append(selected_node)
                Network[selected_node].append(new_node)

        # Update the auxiliary list
        aux_list.extend([new_node] * m)  # Add the new node m times
        for node in connected_nodes:
            aux_list.append(node)  # Add each connected node once

    return Network

# Example usage
N = 5  # Total number of nodes
m = 3   # Minimum degree (number of edges per new node)
network = BA(N, m)
print(network)

import matplotlib.pyplot as plt
import numpy as np
from collections import Counter


def generate_BA(N,m):
    
    # Initialize Library
    Network = {node: [] for node in range(N)}



    # Initial nucleus size
    n0 = m + 1

    # Connect every node in the nucleus to each other
    for node1 in range(n0):
            for node2 in range(node1 + 1 , n0):  # More efficient than adding (if node1 != node2) loop
                Network[node1].append(node2)
                Network[node2].append(node1)

    # Debug 
    # print(Network)

    # Initialize Aux. list
    aux_list = []

    # Initialize initial nucleus nodes
    for node in range(n0):
        for i in range(len(Network[node])):
            aux_list.append(node)

    # print(aux_list)


    # Add new nodes from n0 to N (where j is the new node)
    for j in range(n0, N):

        connected_nodes = []
        while len(connected_nodes) < m:
            # Randomly select attachment node
            the_chosen_node = np.random.choice(aux_list)
            if the_chosen_node not in connected_nodes:
                connected_nodes.append(the_chosen_node)

                # Append 'the_chosen_node' at the index 'new_node' in Network
                # Connect the two using the same method as the initial nucleus
                Network[j].append(the_chosen_node)
                Network[the_chosen_node].append(j)


        # Add number i to each connection made
        for i in connected_nodes:
            aux_list.append(i)

        # Add the new node, j * m for each new vertex i
        aux_list.extend([j] * m)  
        
    return Network


def plot_degree_distribution(Network, m):

    k_list = []
    mf_approx_list = []
    for node in Network:
        k = len(Network[node])
        k_list.append(k)

        # Mean Field Approx
        mf_approx = (2*(m^2)) / (k^3)
        mf_approx_list.append(mf_approx)



    # Step 1: Plot the degree distribution as a histogram

    plt.figure(figsize=(8, 6))
    plt.hist(k_list, bins=range(min(k_list), max(k_list) + 2), edgecolor='black', alpha=0.7, label='Degree Distribution')
    plt.xlabel('Degree (k)')
    plt.ylabel('Frequency')
    plt.title('Degree Frequency Distribution')
    plt.legend()
    plt.show()

    return k_list

def plot_equations(k_list, N, m):
    # Incorporate Eq. (3)

    # Create list of degrees from min to max degree in k_list
    degree = range(min(k_list), max(k_list) + 1)

    # Compute P_inf for each degree k in degree
    P_inf = []
    P_N = []
    for k in degree:
        # print(k)
        P_inf_i = (2 * m * (m + 1)) / (k * (k+1) * (k+2))
        P_inf.append(P_inf_i)

        # Add P_N
        P_N_i = P_inf_i * (k / np.sqrt(N))
        P_N.append(P_N_i)


    # # Plot Probability Distribution P_inf
    # plt.figure()
    # plt.plot(degree, P_inf, marker='o', label="P_inf vs k")
    # plt.xlabel("k")
    # plt.ylabel("P_inf")
    # plt.title("Plot of P_inf vs k")
    # plt.grid(True)
    # plt.legend()
    # plt.show()

    # # Plot Probability Distribution P_N
    # plt.plot(degree, P_N, marker='o', label="P_inf vs k")
    # plt.xlabel("k")
    # plt.ylabel("P_N")
    # plt.title("Plot of P_N vs k")
    # plt.grid(True)
    # plt.legend()
    # plt.show()

    # Plot Probability with Frequency Distribution (Normalized)
    # P_inf
    plt.figure()
    plt.plot(degree, P_inf, marker='o', label="P_inf vs k")
    plt.xlabel("k")
    plt.ylabel("P_inf")
    plt.grid(True)
    plt.legend()

    # P_N
    plt.plot(degree, P_N, marker='o', label="P_inf vs k")
    plt.xlabel("k")
    plt.ylabel("P_inf")
    plt.title("Plot of P_inf vs k")
    plt.grid(True)
    plt.legend()

    # Frequency Distribution
    plt.hist(k_list, bins=range(min(k_list), max(k_list) + 2), density=True, edgecolor='black', alpha=0.7, label='Degree Distribution')
    plt.xlabel('Degree (k)')
    plt.ylabel('Frequency / Approximation')
    plt.title('Degree Distribution with Infinite Network Probability')
    plt.legend()
    plt.show()

    return None


m = 4
N = 100


Network = generate_BA(N, m)
k_list = plot_degree_distribution(Network, m)
plot_equations(k_list, N, m)



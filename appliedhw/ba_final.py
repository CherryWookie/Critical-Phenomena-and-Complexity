
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import random 
from scipy.optimize import curve_fit

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

    # Debug
    # print(aux_list)


    # Add new nodes from n0 to N (where j is the new node)
    for j in range(n0, N):

        connected_nodes = []
        while len(connected_nodes) < m:
            # Randomly select attachment node
            the_chosen_node = random.choice(aux_list)
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

def create_degree_list(Network, m):
    k_list = []
    mf_approx_list = []

    for node in Network:
        k = len(Network[node])  # Degree of the node
        k_list.append(k)

        # Correct Mean Field Approximation
        mf_approx = (2 * (m**2)) / (k**3) if k > 0 else 0  # Avoid division by zero
        mf_approx_list.append(mf_approx)

    return k_list

def plot_probabilities(k_list, N, m):
    # Incorporate Eq. (3) and (4)

    # To handle floats in arrays
    k_list = np.array(k_list)
    degree = np.linspace(k_list.min(), k_list.max(), num=100)  

    # Compute P_inf for each degree k in degree
    P_inf = []
    P_mf = []
    P_N = []
    for k in degree:
    
        # Compute P_inf (Equation 3)
        P_inf_i = (2 * m * (m + 1)) / (k * (k + 1) * (k + 2))
        P_inf.append(P_inf_i)

        # Compute P_mf (Equation 2)
        P_mf_i = (2 * (m**2)) / (k**3) if k > 0 else 0
        P_mf.append(P_mf_i)

        # # Compute P_N (Equation 4)
        # P_N_i = P_inf_i * (k / np.sqrt(N))
        # P_N.append(P_N_i)
    
    # P_inf
    plt.figure()
    plt.plot(degree, P_inf, label="P_inf(k)", linewidth=3, color='skyblue')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("k")
    plt.ylabel("P_inf")
    plt.grid(True)
    plt.legend()

    # P_N
    plt.plot(degree, P_mf, label="P_mf(k)", linewidth=3, color='blue')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("k")
    plt.ylabel("P_mf")
    plt.grid(True)
    plt.legend()

    return P_inf

def plot_averaged_probabilities(hist_avg, bins, N, m):
    
    # Degree array
    degree = np.array(bins[:-1])
    
    # Find the first non-zero bin in the hist_avg (Un-comment to optimize graph limits)
    first_non_zero_index = np.argmax(hist_avg > 0)

    # Find first non-zero degree to start Probability calculation at
    trimmed_degree = degree[first_non_zero_index:]
    # trimmed_hist_avg = hist_avg[first_non_zero_index:]

    # Compute trimmed P_inf and P_N
    P_inf = []
    P_mf = []
    P_N = []
    
    for k in trimmed_degree:
        if k > 0:  # Avoid division by zero for k = 0
            # Compute P_inf (Equation 3)
            P_inf_i = (2 * m * (m + 1)) / (k * (k + 1) * (k + 2))
            P_inf.append(P_inf_i)

            # Compute P_mf (Equation 2)
            P_mf_i = (2 * (m**2)) / (k**3) if k > 0 else 0
            P_mf.append(P_mf_i)

            # # Compute P_N (Equation 4)
            # P_N_i = P_inf_i * (k / np.sqrt(N))
            # P_N.append(P_N_i)
        else:
            P_inf.append(0)
            P_N.append(0)

    
    # Plot P_inf and P_N
    plt.figure()
    plt.plot(trimmed_degree, P_inf, label="P_inf(k)", linewidth=3, color='skyblue')
    plt.plot(trimmed_degree, P_mf, label="P_mf(k)", linewidth=3, color='blue')
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel("k (Degree)")
    plt.ylabel("Probability")
    plt.legend()
    plt.grid(True)
    
    ## Optional Axis Limits
    # plt.xlim(trimmed_degree[0], trimmed_degree[-1])

def compute_p_inf(hist_avg, bins, m):

    degree = np.array(bins[:-1])

    # Compute Pinfvals
    P_inf_vals = []
    for k in degree:
        if k > 0:  # Avoid division by zero for k = 0
            # Compute P_inf (Equation 3)
            P_inf_i = (2 * m * (m + 1)) / (k * (k + 1) * (k + 2))
            P_inf_vals.append(P_inf_i)
        else:
            P_inf_vals.append(0)
            
    return np.array(P_inf_vals)

def plot_degree_distribution(k_list, N, m):


    # Convert to array
    k_list = np.array(k_list)

    # Compute unique degrees frequencies
    unique_degrees, counts = np.unique(k_list, return_counts=True)

    # Normalize frequencies to obtain probabilities
    probabilities = counts / len(k_list)

    # Plot the degree distribution as discrete points
    plt.plot(unique_degrees, probabilities, 'o', markersize=3, label='Degree Distribution', color='black')
    
    # Add labels, title, and legend
    plt.xlabel('Degree (k)')
    plt.ylabel('Probability')
    plt.xscale('log')
    plt.yscale('log')
    plt.title(f'Degree Distribution with Infinite Network Probability\nNetwork Size: {N}')
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_averaged_degree_distribution(hist_avg, N, m):

    hist_avg = np.array(hist_avg)

    degrees = np.arange(len(hist_avg))

    hist_avg_normalized = hist_avg / np.sum(hist_avg)

    lowest_value = degrees[np.nonzero(hist_avg_normalized)[0][0]]


    # Plot averaged degree distribution as discrete points
    plt.plot(degrees, hist_avg_normalized, 'o', markersize=4, label='Averaged Degree Distribution', color='black')
    plt.xlim(left=lowest_value - 1)
    plt.xlabel('Degree (k)')
    plt.ylabel('Probability')
    plt.xscale('log')
    plt.yscale('log')
    plt.title(f'Averaged Degree Distribution with Infinite Network Probability\nNetwork Size: {N}')
    plt.legend()
    plt.grid(True)
    plt.show()

# %%
def plot_finite_size_correction(hist_avg, bins, N_vals, m, M):
    plt.figure()
    for N in N_vals:
        # Initialize sum of histograms for averaging
        hist_sum = None

        # Initialize Max Degree Value
        max_degree = 0
    
        for network in range(M):
            Network_i = generate_BA(N, m)
            k_list = create_degree_list(Network_i, m)

            # Keep track of max_degree encountered
            max_degree = max(max_degree, max(k_list))

        # Define consitent bins for np.histogram()
        bins = range(0, max_degree + 2)

        for network in range(M):
            Network_i = generate_BA(N, m)
            k_list = create_degree_list(Network_i, m)

            hist, _ = np.histogram(k_list, bins=bins) # Need to use consistent bins here

            if hist_sum is None:
                hist_sum = hist
            else:
                hist_sum += hist
                
        # Normalization
        hist_avg = hist_sum / M
        hist_avg = hist_avg / np.sum(hist_avg)

        # Degree
        degree = np.array(bins[:-1])

        # Compute the theoretical P_inf values for this network size
        P_inf_vals = compute_p_inf(hist_avg,bins,m) # Outputs an array of P_inf (theoretical values)

        correction = hist_avg / P_inf_vals  # Finite size correction

        # Plot the correction
        plt.plot(degree, correction, label=f'N = {N}', linewidth=2)

        # Enable grid lines
        plt.grid(which='major', linestyle='-', linewidth=0.5, color='black')  # Major grid
        plt.grid(which='minor', linestyle='--', linewidth=0.5, color='gray')  # Minor grid

        # Turn on minor ticks
        plt.minorticks_on()

        # Labels
        plt.xlabel('Degree k')
        plt.ylabel(r'$P_{\text{est}}(k) / P_{\infty}(k)$')
        plt.xscale('log')
        plt.yscale('log')
        # plt.grid(True)
        plt.legend()
        plt.title('Finite Size Corrections: Pest(k) / P∞(k)')

    plt.show()

# Repeat most of the previous function but solving for w(x)
def plot_Wx(hist_avg, bins, N_vals, m, M):
    plt.figure()
    all_x = [] # relationship x = k/sqrt(N)
    all_w = [] # function w(x)

    for N in N_vals:
        # Initialize sum of histograms for averaging
        hist_sum = None

        # Initialize Max Degree Value
        max_degree = 0
    
        for network in range(M):
            Network_i = generate_BA(N, m)
            k_list = create_degree_list(Network_i, m)

            # Keep track of max_degree encountered
            max_degree = max(max_degree, max(k_list))

        # Define consitent bins for np.histogram()
        bins = range(0, max_degree + 2)

        for network in range(M):
            Network_i = generate_BA(N, m)
            k_list = create_degree_list(Network_i, m)

            hist, _ = np.histogram(k_list, bins=bins) # Need to use consistent bins here

            if hist_sum is None:
                hist_sum = hist
            else:
                hist_sum += hist
                
        # Normalization
        hist_avg = hist_sum / M
        hist_avg = hist_avg / np.sum(hist_avg)

        # Degree
        degree = np.array(bins[:-1])

        # Compute the theoretical P_inf values for this network size
        P_inf_vals = compute_p_inf(hist_avg,bins,m) # Outputs an array of P_inf (theoretical values)

        correction = hist_avg / P_inf_vals  # Finite size correction

        # Compute x
        x = degree / np.sqrt(N)
        all_x.extend(x)
        all_w.extend(correction)

    all_x = np.array(all_x)
    all_w = np.array(all_w)

    # Plot w(x)
    plt.scatter(all_x, all_w, alpha=0.4, label='w(x)')


    # Enable grid lines
    plt.grid(which='major', linestyle='-', linewidth=0.5, color='black')  # Major grid
    plt.grid(which='minor', linestyle='--', linewidth=0.5, color='gray')  # Minor grid

    # Turn on minor ticks
    plt.minorticks_on()

    # Labels
    plt.xlabel(r'$x = \frac{k}{\sqrt{N}}$')
    plt.ylabel(r'$w(x)$')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.title('Universal Function w(x)')
    plt.legend()

    plt.show()

def compute_avg_cluster_coef(N_vals,M,m):

    # Compute average cluster coefficient
    avg_cc = []

    for N in N_vals:
        for _ in range(M):
            Network_i= generate_BA(N,m)
            k_list = create_degree_list(Network_i, m)
            cc_sum = 0
            for node in Network_i:
                k = k_list[node]
                connections = Network_i[node] # connections between the node itself and others
                if k < 2:
                    c_i = 0 # k has to be bigger than 2 to form a connection
                else:
                    e_i = 0
                    for i in connections:
                        for j in connections:
                            if i in Network_i[j]:
                                e_i += 1   # Plus one edge
                    c_i = (2 * e_i) / (k * (k - 1))
                cc_sum += c_i
            
            avg_cc_i = (1 / N) * cc_sum
        avg_cc.append(avg_cc_i / M)

    # Inverse N plot for comparison
    x = np.linspace(min(N_vals), max(N_vals), 500 )
    y = 1 / x

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(N_vals, avg_cc, marker='o', linestyle='-', color='r', label="Average Clustering Coefficient")
    plt.plot(x, y, linestyle='-', color='gray', label='1/N Power Law')
    plt.xlabel("Network Size (N)", fontsize=12)
    plt.ylabel("Average Clustering Coefficient", fontsize=12)
    plt.title("Average Clustering Coefficient vs Network Size (N)", fontsize=14)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # Log-Log Plot
    plt.figure(figsize=(8, 6))
    plt.plot(N_vals, avg_cc, marker='o', linestyle='-', color='r', label="Average Clustering Coefficient")
    plt.plot(x, y, linestyle='-', color='gray', label='1/N Power Law')
    plt.xlabel("Network Size (N)", fontsize=12)
    plt.ylabel("Average Clustering Coefficient", fontsize=12)
    plt.title("Average Clustering Coefficient vs Network Size (N) LOG-LOG Scale", fontsize=14)
    plt.grid(True)
    plt.xscale('log')  
    plt.yscale('log')  
    plt.legend()
    plt.tight_layout()
    plt.show()

    return avg_cc

m = 4
N = 10**6

np.random.seed(9)

# Compute Network for N = 10^6 nodes
Network = generate_BA(N, m)
k_list = create_degree_list(Network, m)
plot_probabilities(k_list, N, m)
plot_degree_distribution(k_list, N, m)

# For M different Networks
M = 10000
N_vals = [50, 100, 200, 500]

for N in N_vals:

    hist_sum = None

    # Initialize Max Degree Value
    max_degree = 0

    for network in range(M):
        Network_i = generate_BA(N, m)
        k_list = create_degree_list(Network_i, m)

        # Keep track of max_degree encountered
        max_degree = max(max_degree, max(k_list))

    # Define consitent bins for np.histogram()
    bins = range(0, max_degree + 2)

    for network in range(M):
        Network_i = generate_BA(N, m)
        k_list = create_degree_list(Network_i, m)

        hist, _ = np.histogram(k_list, bins=bins) # Need to use consistent bins here

        if hist_sum is None:
            hist_sum = hist
        else:
            hist_sum += hist
            
    # Normalization
    hist_avg = hist_sum / M
    hist_avg = hist_avg / np.sum(hist_avg)

    # Plot the histograms
    plot_averaged_probabilities(hist_avg, bins, N, m)
    plot_averaged_degree_distribution(hist_avg, N, m)


# Plot the Finite Size Correction   
plot_finite_size_correction(hist_avg, bins, N_vals, m, M)
plot_Wx(hist_avg, bins, N_vals, m, M)

m = 4
N_vals = [25, 100, 200, 500, 1000, 5000, 10000]
M = 25

compute_avg_cluster_coef(N_vals,M,m)



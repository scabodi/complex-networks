from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import random as rn

# ====================== FUNCTIONS USED BY THE MAIN CODE ===================
#
# ====================== FOR THE MAIN CODE SCROLL TO THE BOTTOM ============

def get_giant_size(net):
    """
    Calculates the size of the largest component (i.e. the giant component) of
    the network.

    Parameters
    ----------
    net: networkx.Graph() object

    Returns
    -------
    giant_size: int
        size of the giant component

    """
    components = nx.connected_components(net)
    giant_size = len(max(components, key=len))
    return giant_size


def simulate_edge_removal(orignet, order):
    """
    Performs an edge removal simulation

    Parameters
    ----------
    orignet: networkx.Graph() object
        Network in which the edge removal is simulated. A copy of orignet is
        created for the simulations, and the original network is not changed.
    order: list of tuples
        network edges sorted in the order in which they will be removed

    Returns
    -------
    giant_sizes: np.array of ints
        sizes of the giant component at different edge densities
    """
    giant_sizes = []
    net = orignet.copy() # Creating a copy of the original network
    n = len(orignet.edges())
    # YOUR CODE HERE
    #TODO: Loop over edges and remove them in given order.
    for edge in order:
        if net.has_edge(*edge):
            net.remove_edge(*edge)
            giant_sizes.append(get_giant_size(net))
    return giant_sizes

def run_link_removal(path, net_name):
    """
    Sets up framework and runs the edge removal simulation.

    Parameters
    ----------
    path: string
        path to the network to be analyzed
    net_name: string
        name of the network (for labeling)

    Returns
    -------
    No direct output, saves figure of the giant component size as a function
    of network density.
    """
    # setting up:
    # YOUR CODE HERE
    net = nx.read_weighted_edgelist(path)
    #net = nx.Graph() # Read the network from path
    N = len(net.nodes()) # Replace with the number of nodes
    edges = net.edges() # Replace with the network edges
    # random = [rn.randint(0, len(edges)) for i in range(len(edges))]
    # nx.set_edge_attributes(net, random, 'rand_num')
    tot_edges = len(edges)
    for e in edges:
        net[e[0]][e[1]]['rand_num'] = rn.randint(0, tot_edges)
        print(net.get_edge_data(*e).get('rand_num'))

    fig = plt.figure(figsize=(16, 16 * 3 / 4.))
    ax = fig.add_subplot(111)
    fig.suptitle(net_name)

    # defining orders in which to remove the edges

    # YOUR CODE HERE
    descending_weight_edge_order = sorted(edges, key=lambda edge:net.get_edge_data(*edge).get('weight'), reverse=True) # edges sorted by decreasing weight
    ascending_weight_edge_order = sorted(edges, key=lambda edge:net.get_edge_data(*edge).get('weight')) # edges sorted by increasing weight
    random_edge_order = sorted(edges, key=lambda edge:net.get_edge_data(*edge).get('rand_num')) # edges sorted in random order

    print('computing betweenness')
    # YOUR CODE HERE
    edge_to_ebc = nx.edge_betweenness_centrality(net) # Replace with a dictionary of edge betweennes values
    print('ended')

    # sorting the edges by their betweenness:
    # YOUR CODE HERE
    sorted_list = sorted(edge_to_ebc.items(), key=lambda x: x[1], reverse=True)
    ebc_edge_order = [x[0] for x in sorted_list]
    #TODO: Replace by edges sorted by decreasing edge betweenness, i.e. sort the dictionary keys by the values

    # edge removal:

    for order, order_name, color, ls, lw in zip(
        [descending_weight_edge_order, ascending_weight_edge_order,
         random_edge_order, ebc_edge_order],
        ["w_big_first",
         "w_small_first", 'random', "betweenness"],
        ["r", "y", "b", "k"],
        ["-", "-", "-", "-"],
        [2, 3, 4, 5]):

        print(order_name)

        giant_sizes = simulate_edge_removal(net, order)
        fracs = np.linspace(0, 1, len(giant_sizes))

        ax.plot(fracs, np.array(giant_sizes)/ float(N), "-", color=color, ls=ls,
                label="g " + order_name, lw=lw)

        # YOUR CODE HERE
        ax.set_ylabel('Size of largest component S') # Set label
        ax.set_xlabel('Fraction of removed links f') # Set label

        ax.legend(loc=2)

    return fig

# =========================== MAIN CODE BELOW ==============================

if __name__ == "__main__":

    network_path = './OClinks_w_undir.edg' # You may want to change the path to the edge list file
    network_name = 'fb-like-network'

    fig = run_link_removal(network_path, network_name)
    fig.savefig("./fb_like_error_and_attack_tolerance.pdf")

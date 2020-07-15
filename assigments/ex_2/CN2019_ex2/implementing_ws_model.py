from __future__ import print_function
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random

# ====================== FUNCTIONS USED BY THE MAIN CODE ===================
#
# ====================== FOR THE MAIN CODE SCROLL TO THE BOTTOM ============
def next(n, i, m):
    return (i+m)%n
def prev(n, i, m):
    return (i+n-m)%n

def ring(n, m):
    """
    This function creates the basic ring (to be rewired) with n nodes
    in which each node is connected to m nodes on the left and right.

    Parameters
    ----------
    n : int
      Number of nodes
    m : int
      Number of neighbors to connect left and right

    Returns
    -------
    network : graph
             The basic ring before rewiring
    """
    network = nx.Graph()

    for i in range(n):
        network.add_node(i)
    for i in range(n):
        edges = []
        for j in range(m):
            edges.append((i, next(n, i, j+1)))
            edges.append((i, prev(n, i, j+1)))
        network.add_edges_from(edges)

    nx.draw_circular(network)
    return network

def ws(n, m, p):
    """
    This function call the ring() function to make a basic ring and then
    rewires each link with  probability p and also prints the total number of
    links and the number of rewired links.
    Note self-loops are not allowed when rewiring (check that you do not rewire
    the end of a link to the node at its other end!)

    Parameters
    ----------
    n : int
      Number of nodes
    m : int
      Number of neighbors to connect left and right
    p : float
        Rewiring probability

    Returns
    -------
    network : graph
        The Watts-Strogatz small-world network

    """
    network = ring(n, m)
    edges = network.edges()
    nodes = network.nodes()
    rewired_num = 0 # tracks the number of rewired links
    total_num = len(network.edges()) # tracks the total number of links in the network

    for e in edges:
        if np.random.rand() < p: #rewire edge
            n = e[0]
            non_neighbors = list(nx.non_neighbors(network, n))
            new_node = np.random.choice(non_neighbors)
            non_neighbors.remove(new_node)
            network.remove_edge(*e)
            network.add_edge(n , new_node)
            rewired_num += 1

    print("total number of links:")
    print(total_num)
    print("number of rewired links:")
    print(rewired_num)
    return network

# =========================== MAIN CODE BELOW ==============================

if __name__ == "__main__":
    np.random.seed(42)
    #visualizing the rings for p = 0 ...
    graph1 = ws(15, 2, .1)
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    nx.draw_circular(graph1)

    figure_filename = 'Small_world_ring.pdf'
    fig1.savefig(figure_filename)

    total_num_edges = len(list(graph1.edges()))
    print("Total number of edges for n = 15, m = 2, p = 0 :")
    print(total_num_edges)
    #... and p = 0.5
    graph2 = ws(100, 2, 0.5)
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    nx.draw_circular(graph2, node_size=20)

    figure_filename = './Small_world_rewired.pdf'

    fig2.savefig(figure_filename)

    # Produce the basic ring network and calculate the average clustering
    # coefficient and average shortest path of the network
    basic_ring = ring(500, 2)

    c_basic = nx.average_clustering(basic_ring, count_zeros=True)
    l_basic = nx.average_shortest_path_length(basic_ring)

    probability = [0.001*(2**n) for n in range(11)] #[0.001, 0.002, 0.004, ...]
    relative_c = []
    relative_l = []

    for p in probability:
        smallworld = ws(500, 2, p)

        # gets all connected components; mostly there is just one:
        components = nx.connected_component_subgraphs(smallworld)
        # finds the largest to be used for the average path length:
        largest_component = max(components, key=len)

        c_rewired = nx.average_clustering(smallworld, count_zeros=True) / c_basic
        l_rewired = nx.average_shortest_path_length(largest_component) / l_basic
        # Update relative_c and relative_l
        relative_c.append(c_rewired)
        relative_l.append(l_rewired)

    fig3 = plt.figure()
    ax = fig3.add_subplot(111)
    ax.semilogx(probability, relative_c, marker='o', ls='-', color='b',
                label='relative avg c')
    ax.semilogx(probability, relative_l, marker='o', ls='-', color='r',
                label='relative avg shortest path')

    # Label the axes
    ax.set_xlabel('p') # sets the label of the x axis
    ax.set_ylabel('Relative quantities')
    ax.legend()
    figure_filename = 'WS_relative_c_and_l.pdf'

    fig3.savefig(figure_filename)

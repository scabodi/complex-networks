from __future__ import print_function
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import binned_statistic

# ====================== FUNCTIONS USED BY THE MAIN CODE ===================
#
# ====================== FOR THE MAIN CODE SCROLL TO THE BOTTOM ============

def lin_log_bins(max_degree):
    # lin-log binning: for k=1..10 use linear bins, then logarithmic bins
    # have the number of logbins such that there are 10 bins per decade

    num_logbins = int(np.log10(1.5 * max_degree) - np.log10(1.5)) * 10

    # generate log bins from k=1.5 to k=1.5*max(degree)
    bins = np.logspace(
        np.log10(1.5), np.log10(1.5 * max_degree), num_logbins)

    return bins


def rand_prob_node(G, m):
    nodes_probs = [] 
    degrees = [val for (node, val) in G.degree()]
    sum_degrees = sum(degrees)
    for node in G.nodes():
        node_degr = G.degree(node)
        node_proba = node_degr / sum_degrees
        nodes_probs.append(node_proba)
    random_proba_node = np.random.choice(G.nodes(), p=nodes_probs, size=m, replace=False)
    return random_proba_node


def ba_network(n, m, seedsize=3):
    # Generate initial small seed network (clique of seedside nodes)
    net = nx.complete_graph(seedsize)
    # YOUR CODE HERE
    # Grow the network here
    nodes = net.nodes()
    i = len(nodes)

    while i < n: #continue growing until I have N nodes
        node = rand_prob_node(net, m)
        net.add_node(i)
        for j in range(m):
            net.add_edge(i, node[j])
        i += 1

    return net

# =========================== MAIN CODE BELOW ==============================

if __name__ == "__main__":
    np.random.seed(42)

    # part a
    fig = plt.figure()
    ax = fig.add_subplot(111)

    net = ba_network(200, 1)
    nodes = net.nodes()
    degrees_dict = nx.degree(net)
    degrees = [degrees_dict[node] for node in nodes]

    print("The maximum degree is: ", max(degrees))
    print("The total number of edges is: ", len(net.edges()))

    nx.draw_spring(
        net, node_size=100, node_color=degrees, cmap='autumn',
        vmin=np.min(degrees), vmax=np.max(degrees))
    ax.set_aspect('equal')

    figure_filename = 'BA_visualized.pdf'

    fig.savefig(figure_filename)
    # or just use plt.show() and save manually

    # part b
    net = ba_network(10000, 2)
    degrees = [deg for _, deg in nx.degree(net)]
    # if you are using an older version of networkx where the return value of nx.degree is a dict instead of
    # a DegreeView, you will get a type error from the above line. To fix, change it to:
    # degrees = list(nx.degree(net).values())

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # so use np.histogram to get histogram and bin edges
    bins = lin_log_bins(max(degrees))
    pk, bin_edges = np.histogram(degrees, bins=bins, density=True)

    bincenters, _, _ = binned_statistic(
        degrees, degrees, statistic='mean', bins=bins)
    ax.set_xlabel('Degree k')
    ax.set_ylabel('P(k)')

    ax.loglog(bincenters, pk, 'ro', label='Simulated')
    ax.loglog(bins, 2 * 2 * (2 + 1) /
              (bins * (bins + 1) * (bins + 2)),
              label='Theoretical')

    ax.legend()

    figure_filename = 'BA_degree_distribution.pdf'

    fig.savefig(figure_filename)
    # or just use plt.show() and save manually

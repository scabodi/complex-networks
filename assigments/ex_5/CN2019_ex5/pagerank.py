# In order to understand how the code works, it is a good idea to check the
# final section of the file that starts with
#   if __name__ == '__main__'
#
# Your task is essentially to replace all the parts marked as TODO or as
# instructed through comments in the functions, as well as the filenames
# and labels in the main part of the code which can be found at the end of
# this file.
#
# The raise command is used to help you out in finding where you still need to
# write your own code. When you successfully modified the code in that part,
# remove the `raise` command.
from __future__ import print_function
import timeit

import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import networkx as nx
import random as rd
import time

# ====================== FUNCTIONS USED BY THE MAIN CODE ===================
#
# ====================== FOR THE MAIN CODE SCROLL TO THE BOTTOM ============

###############################################################
# Code that is given to you, and does not need to be modified #
###############################################################

def add_colorbar(cvalues, cmap='OrRd', cb_ax=None):
    """
    Add a colorbar to the axes.

    Parameters
    ----------
    cvalues : 1D array of floats

    """
    eps = np.maximum(0.0000000001, np.min(cvalues)/1000.)
    vmin = np.min(cvalues) - eps
    vmax = np.max(cvalues)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    scm = mpl.cm.ScalarMappable(norm, cmap)
    scm.set_array(cvalues)
    if cb_ax is None:
        plt.colorbar(scm)
    else:
        cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, orientation='vertical')

######################################################
# Starting from here you might need to edit the code #
######################################################

def pageRank(network, d, n_steps):
    """
    Returns the PageRank value, calculated using a random walker, of each
    node in the network. The random walker teleports to a random node with
    probability 1-d and with probability d picks one of the neighbors of the
    current node.

    Parameters
    -----------
    network : a networkx graph object
    d : damping factor of the simulation
    n_steps : number of steps of the random walker

    Returns
    --------
    page_rank: dictionary of node PageRank values (fraction of time spent in
               each node)
    """

    # Initializing PageRank dictionary:
    pageRank = {}
    nodes = list(network.nodes())

    # YOUR CODE HERE
    #TODO: write code for calculating PageRank of each node
    # 1) Initialize PageRank of each node to 0
    # 2) Pick a random starting point for the random walker (Hint: np.random.choice)
    # 3) Random walker steps, at each step:
    #   1) Increase the PageRank of current node by 1
    #   2) Check if the random walker will teleport or go to a neighbor
    #   3) Pick the next node either randomly or from the neighbors
    #   4) Update the current node variable to point to the next node
    # 4) Repeat random walker steps 1-4 n_steps times
    # 5) Normalize PageRank by n_steps

    for node in nodes:
      pageRank[node] = 0

    startNode = np.random.choice(nodes)
    currNode = startNode

    for step in range(n_steps):
        pageRank[currNode] += 1
        r = rd.random()
        if (r < d):
            neighbors = list(nx.neighbors(network,currNode))
            currNode = np.random.choice(neighbors)
        else:
            currNode = np.random.choice(nodes)

    for node in nodes:
      pageRank[node] /= n_steps

    return pageRank

def pagerank_poweriter(g, d, iterations):
    """
    Uses the power iteration method to calculate PageRank value for each node
    in the network.

    Parameters
    -----------
    g : a networkx graph object
    d : damping factor of the simulation
    iterations : number of iterations to perform

    Returns
    --------
    pr_old : dict where keys are nodes and values are PageRank values
    """
    print("Running function for obtaining PageRank by power iteration...")
    # YOUR CODE HERE
    #TODO: write code for calculating power iteration PageRank
    # Some pseudocode:
    # 1) Create a PageRank dictionary and initialize the PageRank of each node
    #    to 1/n where n is the number of nodes.
    # 2) For each node i, find nodes having directed link to i and calculate
    #    sum(x_j(t-1)/k_j^out) where the sum is across the node's neighbors
    #    and x_j(t-1) is the PageRank of node  .
    # 3) Update each node's PageRank to (1-d)*1/n + d*sum(x_j(t-1)/k_j^out).
    # 4) Repeat 2-3 n_iterations times.

    pr_old, pr_new = {}, {}
    nodes = list(network.nodes())
    n = len(nodes)

    initialRank = 1/n
    for node in nodes:
        pr_old[node] = initialRank
        pr_new[node] = initialRank

    for t in range(1, iterations):
        for i in nodes:
            incoming_edges = g.in_edges(i)
            v_i = [e[0] for e in incoming_edges]  #take starting node for incoming edges of node i
            sum = 0
            for j in v_i:
                k_j = len(g.out_edges(j))
                sum += pr_old[j]/k_j
            pr_new[i] = (1-d)*(1/n) + d*sum
        # sanity checks in each itteration:
        '''print('PageRank sums to ')
        print(sum(pr_new.values()))
        print('PageRank difference since last iteration:')
        print(sum([abs(pr_new[i]-pr_old[i]) for i in g]))'''
        for i in nodes:
            pr_old[i] = pr_new[i]
    return pr_old

def visualize_network(network, node_positions, cmap='OrRd',
                      node_size=3000, node_colors=[], with_labels=True, title=""):
    """
    Visualizes the given network using networkx.draw and saves it to the given
    path.

    Parameters
    ----------
    network : a networkx graph object
    node_positions : a dict of positions of nodes, obtained by e.g. networkx.graphviz_layout
    cmap : colormap
    node_size : int
    node_colors : a list of node colors
    with_labels : should node labels be drawn or not, boolean
    title: title of the figure, string
    """
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)
    if node_colors:
        #visualize the networks with nodes colored by PageRank.
        nx.draw(network, pos=node_positions, node_color = node_colors, cmap=cmap, node_size=node_size, with_labels=with_labels)
        add_colorbar(node_colors)
    else:
        #visualize the networks without node coloring.8
        nx.draw(network, pos=node_positions, cmap=cmap, node_size=node_size, with_labels=with_labels)
    ax.set_title(title)
    plt.tight_layout()

    return fig

def investigate_d(network, ds, colors, n_steps):
    """
    Calculates PageRank at different values of the damping factor d and
    visualizes and saves results for interpretation

    Parameters
    ----------
    network : a NetworkX graph object
    ds : a list of d values
    colors : visualization color for PageRank at each d, must have same length as ds
    n_steps : int; number of steps taken in random walker algorithm
    """
    #import pdb; pdb.set_trace()
    nodes = network.nodes()
    n_nodes = len(nodes)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #TODO: write a for loop to obtain node PageRank values at each d and to plot the PageRank.
    # YOUR CODE HERE
    # use zip to loop over ds and colors at once
    ds_colors = zip(ds, colors)
    for d,c in ds_colors:
        d = round(d, 1)
        pageRank_rw = pageRank(network, d, n_steps)
        pageRank_rw_array = np.zeros(n_nodes)
        for node in nodes:
            pageRank_rw_array[int(node)] = pageRank_rw[node]
        plt.plot(range(0, n_nodes), pageRank_rw_array, color=c, label=r'Rw d='+str(d))

    #plt.plot(range(0, n_nodes), mapped, 'k+', label=r'Random walker')
    ax.set_xlabel(r'Node index')
    ax.set_ylabel(r'PageRank')
    ax.set_title(r'PageRank with different damping factors')
    ax.legend(loc=0)
    plt.tight_layout

    return fig

# =========================== MAIN CODE BELOW ==============================

if __name__ == '__main__':

    #TODO: replace, set the correct path to the network
    network_path = './pagerank_network.edg'

    network =  nx.read_edgelist(network_path, create_using=nx.DiGraph())

    # Visualization of the network (note that spring_layout
    # is intended to be used with undirected networks):
    node_positions = nx.spring_layout(network.to_undirected())
    cmap = 'OrRd'
    node_size = 3000

    fig=visualize_network(network, node_positions, cmap=cmap, node_size=node_size, title="Network")

    fig.savefig('./network_visualization.pdf')

    nodes = network.nodes()
    n_nodes = len(nodes)
    # PageRank with self-written function
    # YOUR CODE HERE
    n_steps = 10000 # TODO: replace: set a reasonable n_steps
    d = 0.85 # TODO: replace: set a reasonable d; nx.pagerank uses d = 0.85
    pageRank_rw = pageRank(network, d, n_steps)

    # Visualization of PageRank on network:
    node_colors = [pageRank_rw[node] for node in nodes]

    fig = visualize_network(network, node_positions, cmap=cmap, node_size=node_size,
                            node_colors=node_colors, title="PageRank random walker")

    fig.savefig('./network_visualization_pagerank_rw.pdf')

    # PageRank with networkx:
    pageRank_nx = nx.pagerank(network)

    # Visualization to check that results from own function and nx.pagerank match:

    pageRank_rw_array = np.zeros(n_nodes)
    pageRank_nx_array = np.zeros(n_nodes)
    for node in nodes:
        pageRank_rw_array[int(node)] = pageRank_rw[node]
        pageRank_nx_array[int(node)] = pageRank_nx[node]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(range(0, n_nodes), pageRank_rw_array, 'k+', label=r'Random walker')
    plt.plot(range(0, n_nodes), pageRank_nx_array, 'rx', label='networkx')
    ax.set_xlabel(r'Node index')
    ax.set_ylabel(r'PageRank')
    ax.set_title(r'PageRank with different methods')
    ax.legend(loc=0)
    plt.tight_layout()

    fig.savefig('./network_visualization_pagerank_nx.pdf')

    # PageRank with power iteration
    n_iterations = 10
    pageRank_pi = pagerank_poweriter(network, d, n_iterations)

    # Visualization of PageRank by power iteration

    node_colors = [pageRank_pi[node] for node in nodes]
    fig = visualize_network(network, node_positions, cmap=cmap, node_size=node_size,
                            node_colors=node_colors, title="PageRank power iteration")

    fig.savefig('./network_visualization_pagerank_pi.pdf')

    # Investigating th  of the power iteration fuction
    num_tests = 3 
    # YOUR CODE HERE
    k5net = nx.directed_configuration_model(10**3*[5], 10**3*[5],create_using = nx.DiGraph)
    # TODO: Print results: how many seconds were taken for the test network of
    # 10**4 nodes, how many hours would a 26*10**6 nodes network take?
    # Run num_tests times and use the average or minimum running time in your calculations.
    t0 = time.clock()
    # Investigating the running time of the random walker function
    n_nodes = 10**3
    # YOUR CODE HERE
    n_steps = 1000*n_nodes # TODO: set such number of steps that each node gets visited on average 1000 times
    pageRank_rw = pageRank(k5net, d, n_steps)
    print('time elapsed for random : ', (time.clock()-t0))

    t0 = time.clock()
    pageRank_pi = pagerank_poweriter(k5net, d, n_nodes*10)
    print('time elapsed for power iteration : ', (time.clock()-t0))

    # Investigating effects of d:
    ds = np.arange(0, 1.2, 0.2)
    colors = ['b', 'r', 'g', 'm', 'k', 'c']

    fig = investigate_d(network, ds, colors, 10000)

    fig.savefig('./investigate_role_of_damping_factor.pdf')

    # YOUR CODE HERE
    network_path_wp = './wikipedia_network.edg' # TODO: replace, set the correct network path
    network_wp = nx.read_edgelist(network_path_wp, create_using=nx.DiGraph())
    # TODO: replace with the network loaded with nx.read_edgelist.

    pageRank_wp = dict(nx.pagerank(network_wp))
    indegree_wp = dict(network_wp.in_degree())
    outdegree_wp = dict(network_wp.out_degree())

    if pageRank_wp is not {}:
        highest_pr = sorted(pageRank_wp, key=lambda k: pageRank_wp[k])[::-1][0:5]
        print('---Highest PageRank:---')
        for p in highest_pr:
            print(pageRank_wp[p], ":", p)
    if indegree_wp is not {}:
        highest_id = sorted(indegree_wp, key=lambda k: indegree_wp[k])[::-1][0:5]
        print('---Highest In-degree:---')
        for p in highest_id:
            print(indegree_wp[p], ":", p)
    if outdegree_wp is not {}:
        highest_od = sorted(outdegree_wp, key=lambda k: outdegree_wp[k])[::-1][0:5]
        print('---Highest Out-degree:---')
        for p in highest_od:
            print(outdegree_wp[p], ":", p)

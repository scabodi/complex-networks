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

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

# ====================== FUNCTIONS USED BY THE MAIN CODE ===================
#
# ====================== FOR THE MAIN CODE SCROLL TO THE BOTTOM ============


###############################################################
# Code that is given to you, and does not need to be modified #
###############################################################

def plot_network_usa(net, xycoords, bg_figname, edges=None, alpha=0.3):
    """
    Plot the network usa.

    Parameters
    ----------
    net : the network to be plotted
    xycoords : dictionary of node_id to coordinates (x,y)
    edges : list of node index tuples (node_i,node_j),
            if None all network edges are plotted.
    alpha : float between 0 and 1, describing the level of
            transparency
    """
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 0.9])
    # ([0, 0, 1, 1])
    img = plt.imread(bg_figname)
    axis_extent = (-6674391.856090588, 4922626.076444283,
                   -2028869.260519173, 4658558.416671531)
    ax.imshow(img, extent=axis_extent)
    ax.set_xlim((axis_extent[0], axis_extent[1]))
    ax.set_ylim((axis_extent[2], axis_extent[3]))
    ax.set_axis_off()
    nx.draw_networkx(net,
                     pos=xycoords,
                     with_labels=False,
                     node_color='k',
                     node_size=5,
                     edge_color='r',
                     alpha=alpha,
                     edgelist=edges)
    return fig, ax
######################################################
# Starting from here you might need to edit the code #
######################################################


# =========================== MAIN CODE BELOW ==============================

if __name__ == '__main__':
    #TIP: #TODO: set correct paths where the initial provided information is
    csv_path = './US_airport_id_info.csv'
    network_path = "./aggregated_US_air_traffic_network_undir.edg"
    bg_figname = './US_air_bg.png'

    id_data = np.genfromtxt(csv_path, delimiter=',', dtype=None, names=True)
    xycoords = {}
    for row in id_data:
        xycoords[str(row['id'])] = (row['xcoordviz'], row['ycoordviz'])
    net = nx.read_weighted_edgelist(network_path)

    #TIP: #TODO: Write your code here to calculate the required metrics, compute the
    #TIP: # thresholded network, maximal and minimal spanning tree, and visualize the networks
    #TIP: # Remember you can use the networkx functions

    # YOUR CODE HERE
    # a) Some network properties:
    # N -- number of network nodes
    N = nx.number_of_nodes(net)
    print('Number of network nodes : '+str(N))
    # L -- number of links
    L = nx.number_of_edges(net)
    print('Number of network links : '+str(L))
    # D -- density
    D = nx.density(net)
    print('Density of the network : '+ str(D))
    # d -- diameter
    max_sub = max(nx.connected_component_subgraphs(net), key=len)
    d = nx.diameter(max_sub)
    print('Diameter of the network : '+str(d))
    # C -- average clustering coefficient
    C = nx.average_clustering(net, count_zeros=True) #TODO verify
    print('Average clustering coefficient of the network : '+str(C))

    # b) plot the network
    fig, ax = plot_network_usa(net, xycoords, bg_figname=bg_figname, edges=net.edges())
    fig.savefig('./US_air_network.png')

    # c) compute min and max spanning tree
    fig1_name = './US_air_max_st.png'
    max_st = nx.maximum_spanning_tree(net)
    fig1, ax = plot_network_usa(net, xycoords, bg_figname=bg_figname, edges=max_st.edges())
    fig1.savefig(fig1_name)

    fig2_name = './US_air_MST.png'
    min_st = nx.minimum_spanning_tree(net)
    fig2, ax = plot_network_usa(net, xycoords, bg_figname=bg_figname, edges=min_st.edges())
    fig2.savefig(fig2_name)

    # d) threashold M = MST number of edges
    fig3_name = './US_air_threshold.png'
    M = nx.number_of_edges(max_st)
    print(str(M))
    #TODO verify that the strongest link are in fact those with weight higher!!
    #order edges based on their weight from highest to lowest
    edges = list(sorted(net.edges(data=True), key=lambda t: t[2].get('weight', 1), reverse=True))
    edges_threshold = edges[:M] # keep those in threshold
    set_th = set([(x[0], x[1]) for x in edges_threshold])
    set_mst = set(max_st.edges())
    print('Number of shared edges = '+str(len(set_mst.intersection(set_th))))
    fig3, ax = plot_network_usa(net, xycoords, bg_figname=bg_figname, edges=edges_threshold)
    fig3.savefig(fig3_name)

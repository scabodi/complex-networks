import sys
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# ====================== FUNCTIONS USED BY THE MAIN CODE ===================
#
# ====================== FOR THE MAIN CODE SCROLL TO THE BOTTOM ============
def get_average_degree(network):
    """
    Returns a list of the degrees of all nodes in the network.

    Parameters
    ----------
    network: a NetworkX graph object

    Returns
    -------
    degrees: list
        degrees of all network node
    """
    degrees = [] # empty list

    N = nx.number_of_nodes(network)
    m = nx.number_of_edges(network)
    nodes = list(network.nodes)


    for n in nodes:
        degrees.append(network.degree(n))
    avg = np.mean(degrees)

    return avg

def ER_properties(n, p):
    '''
    This function builds 100 ER networks with the given parameters and averages
    over them to estimate the expected values of average clustering coefficient,
    average degree and diameter of an ER network with the given parameters. The
    diameter is always computed for the largest ("giant") connected component.

    Parameters
    ----------
    n : int
      Number of nodes
    p : float
      the probability that a pair of nodes are linked is p.

    Returns
    -------
    expected_c: float
                expected value of average clustering coefficient
    expected_k: float
                expected value of average degree
    expected_d: float
    expected value of d*



    Hints:
        The following functions might be useful:
        nx.fast_gnp_random_graph(n, p),
        nx.average_clustering(graph, count_zeros=True),
        nx.connected_component_subgraphs(graph),
        nx.diameter(giant)

        nx.connected_component_subgraphs gives you a list of subgraphs
        To pick the largest, the fastest way is to use max with a key: max(x,key=len)
        returns the longest/largest element of the list x

        For computing averages over realizations, you can e.g. collect your
        values to three lists, c,k,d, and use np.mean to get the average.
    '''

    # YOUR CODE HERE
    n_networks = 100
    c, k, d = [], [], []

    for i in range(n_networks):
        graph = nx.fast_gnp_random_graph(n, p) #create the network with n nodes and p probability
        c.append(nx.average_clustering(graph, count_zeros=True) )#compute the average clustering coefficient for the graph
        k.append(get_average_degree(graph))
        max_sub = max(nx.connected_component_subgraphs(graph), key=len)
        d.append(nx.diameter(max_sub))

    #find the avg value for each list
    expected_c = np.mean(c)
    expected_k = np.mean(k)
    expected_d = np.mean(d)

    return expected_c, expected_k, expected_d

def ER_properties_theoretical(p):
    '''
    This function calculates the theoretical values for clustering coefficients,
    average degree, and diameter for ER networks of size 3 and link probability p.
    The theoretical values can be viewed as expectations, or ensemble averages.
    Therefore, e.g., the expected diameter doesn't have to be integer, although it of
    course always is for a single ER network.

    Parameters
    ----------
    p : float
      the probability that a pair of nodes are linked is p.

    Returns
    -------
    c_theory: float
                theoretical value of average clustering coefficient
    k_theory: float
                theoretical value of average degree
    d_theory: float
                Theoretical value of diameter
    '''

    # Calculate the theoretical values for an ER network with parameters n, p
    n = 3
    c_theory = pow(p, 3)
    k_theory = (n-1)*p
    d_theory = 3*p-2*pow(p, 3)
    return c_theory, k_theory, d_theory


def plot_er_values(n, p_list):
    '''
    This function calculates the theoretical clustering coefficient, average
    degree and diameter for ER network with parameters n and p and plots them
    against the expected values from an ensemble of 100 realizations.

    Parameters
    ----------
    n : int
      Number of nodes
    p_list : list of floats
      where each member is the probability that a pair of nodes are linked.

    Returns
    -------
    fig: matplotlib Figure
                plots of expected values against theoretical values
    '''

    k_list = [] # list for degrees
    c_list = [] # list for clustering coeff values
    d_list = [] # list for diameter values
    c_list_theory = []
    k_list_theory = []
    d_list_theory = []

    for p in p_list:
        print("calculating for n=%d p=%f" % (n, p), file=sys.stderr)
        average_c, average_k, average_d = ER_properties(n, p)
        k_list.append(average_k)
        c_list.append(average_c)
        d_list.append(average_d)

        if n is 3:
            c_theory, k_theory, d_theory = ER_properties_theoretical(p)
            c_list_theory.append(c_theory)
            k_list_theory.append(k_theory)
            d_list_theory.append(d_theory)

    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1, 3, 1) # Three subplots for <K>, <C> and <d*> [(1,3,1) means 1 row, three columns, first subplot)

    ax.plot(p_list, c_list, 'r-', label="Simulated", marker='o')
    if len(c_list_theory) > 0 and c_list_theory[0] is not None:
        ax.plot(p_list, c_list_theory, 'b-', label="Theoretical", marker='o')
    ax.set_xlabel('p')
    ax.set_ylabel('Clustering coefficient')
    ax.legend(loc=0)
    ax.set_title('N = %d' % n, size=18)

    ax2 = fig.add_subplot(1, 3, 2)
    ax2.plot(p_list, k_list, 'r-', label="Simulated", marker='o')
    if len(k_list_theory) > 0 and k_list_theory[0] is not None:
        ax2.plot(p_list, k_list_theory, 'b-', label="Theoretical", marker='o')
    ax2.set_xlabel('p')
    ax2.set_ylabel('Average degree')
    ax2.legend(loc=0)

    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(p_list, d_list, 'r-', label="Simulated", marker='o')
    if len(d_list_theory) > 0 and d_list_theory[0] is not None:
        ax3.plot(p_list, d_list_theory, 'b-', label="Theoretical", marker='o')
    ax3.set_xlabel('p')
    ax3.set_ylabel('Diameter')
    ax3.legend(loc=0)

    fig.tight_layout()

    return fig

# =========================== MAIN CODE BELOW ==============================

if __name__ == "__main__":
    ps = np.arange(0, 1.01, 0.05)

    fig = plot_er_values(n=3, p_list=ps)

    figure_filename = 'ER_properties_3_nodes.pdf'
    fig.savefig(figure_filename)

    # Do the same steps (calculation and visulization) for n=100
    fig = plot_er_values(n=100, p_list=ps)

    figure_filename = 'ER_properties_100_nodes.pdf'
    fig.savefig(figure_filename)

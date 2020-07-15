import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import random as rd
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings("ignore", category=UserWarning)


# ===================== PLOTTING FUNCTIONS ============================

def plot_network_usa(net, xycoords, bg_figname, edges=None, linewidths=None):
    """
    Plot the network usa.

    :param net: the network to be plotted
    :param xycoords: dictionary of node_id to coordinates (x,y)
    :param bg_figname: name of the figure
    :param edges: list of node index tuples (node_i,node_j), if None all network edges are plotted.
    :param linewidths : a list with equal length and order to egdes -list.
    :return: fig: figure, ax: axes
    """

    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 0.9])
    # ([0, 0, 1, 1])
    img = plt.imread(bg_figname)
    axis_extent = (-6674391.856090588, 4922626.076444283, -2028869.260519173, 4658558.416671531)
    ax.imshow(img, extent=axis_extent)
    ax.set_xlim((axis_extent[0], axis_extent[1]))
    ax.set_ylim((axis_extent[2], axis_extent[3]))
    ax.set_axis_off()
    # nx.draw_networkx(net, pos=xycoords, with_labels=False, node_color='k', node_size=5, edge_color='r', alpha=alpha,
    #                 edgelist=edges)
    nx.draw_networkx_nodes(net, pos=xycoords, with_labels=False, node_color='k', node_size=5, alpha=0.2)
    if linewidths is None:
        linewidths = np.ones(len(edges))

    for edge, lw in zip(edges, linewidths):
        nx.draw_networkx_edges(net, pos=xycoords, with_labels=True, edge_color='r', width=lw, edgelist=[edge], alpha=lw)

    return fig, ax


def plot_averaged_prevalence(datavecs, labels, xlabel, ylabel, name, max_x, n_step, colors):
    """
    Plots in a single figure the averaged prevalence of the disease as a function of time for all infection probabilities.

    :param datavecs: data vectors to plot, a list of iterables
    :param labels: labels for the data vectors, list of strings
    :param xlabel: x label for the figure, string
    :param ylabel: y label for the figure, string
    :param name: name of the figure
    :param max_x: max x-value to be plotted
    :param n_step: number of steps used in linspace
    :param colors: list of colors

    :return: fig: figure
    """

    fig = plt.figure(name)
    ax = fig.add_subplot(111)
    x = np.linspace(0, max_x, num=n_step)
    for i in range(len(colors)):
        ax.plot(x, datavecs[i], label=labels[i], color=colors[i], linestyle='-')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc=0)
    ax.grid()

    return fig


def create_scatter_plots(names, x_values, y_values, x_labels, y_label, colors):
    """
    Creates a scatter plot of y_values as a function of x_values.

    :param names: names of the figures
    :param x_values: list of np.arrays
    :param y_values: np.array
    :param x_labels: list of strings
    :param y_label: string
    :param colors: list of strings
    :return: fig: figure
    """
    figures = []

    for name, x_val, x_label, color in zip(names, x_values, x_labels, colors):
        fig = plt.figure(name)
        ax = fig.add_subplot(111)
        ax.plot(x_val, y_values, ls='', color=color, marker='o')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid()
        figures.append(fig)
        plt.close(fig)

    return figures


# ============ SIMULATION AND OTHER FUNCTIONS =============================

def simulate_infection(seed_id, p, immunized=None, infection_edges=None, infection_duration=-1, model='SI'):
    """
    Simulates the infection of the airports

    :param seed_id: int seed node id
    :param p: probability of  a node to be infected
    :param immunized: set of imminized nodes - int values - default value = empty set
    :param infection_edges: dict with key = (airport_a, airport_b) and value = # of times infection passes through link

    for BONUS TASK
    :param model: type of model implemented
    :param infection_duration: duration of the infection - number of seconds

    :return: timeline_number_of_infected: list of the number of infected nodes as function of time
             infection_times: dict with key=node, value=time_of_infection
                             (for the seed takes the first flight's departure time)
    """

    if immunized is None:
        immunized = set()
    if infection_edges is not None:
        infection_source_edge = {}

    timeline_number_of_infected, infection_times = [], {}
    next_infection_times = {}

    # insert the seed node into the dictionary
    infection_times[seed_id] = event_data[0]['StartTime']

    # iterate over all flights
    for row in event_data:
        source = row['Source']
        dest = row['Destination']
        # create the tuple (source, dest) or (dest, source) depending on their ID's order
        if source < dest:
            edge = (source, dest)
        else:
            edge = (dest, source)
        # check if the source node is infected --> is it present in the dictionary?
        source_infected = source in infection_times and row['StartTime'] >= infection_times[source]
        if model == 'SIR' and source_infected and infection_duration > 0:
            if row['StartTime'] >= infection_times[source] + infection_duration:
                # disable infection as  infection_duration elapsed
                source_infected = False
        if model == 'SIS' and source_infected and infection_duration > 0:
            if row['StartTime'] >= infection_times[source] + infection_duration:
                if source not in next_infection_times or row['StartTime'] >= next_infection_times[source] + infection_duration:
                    # disable infection
                    source_infected = False
        if source_infected:
            # check if the target node is susceptible
            dest_susceptible = dest not in infection_times
            if model == 'SIS' and dest in infection_times:
                if row['EndTime'] >= infection_times[dest] + infection_duration:
                    # destinayion healed - susceptible again
                    dest_susceptible = True
            if dest in next_infection_times and row['EndTime'] < next_infection_times[dest] + infection_duration:
                # destinayion healed - susceptible again
                dest_susceptible = False
            if dest_susceptible:
                # infect the target node with probability p
                if rd.random() <= p and (dest not in immunized):
                    if dest not in infection_times:
                        infection_times[dest] = row['EndTime']
                    else:
                        next_infection_times[dest] = row['EndTime']
                    if infection_edges is not None:
                        infection_source_edge[dest] = edge
            else:
                # target node is already infected
                # check if the current flight’s arrival time is earlier than the node’s
                # infection time in the dictionary, and if, update it to the arrival tie
                if row['EndTime'] < infection_times[dest]:
                    infection_times[dest] = row['EndTime']
                    if infection_edges is not None:
                        infection_source_edge[dest] = edge
                else:
                    if dest in next_infection_times and row['EndTime'] < next_infection_times[dest]:
                        next_infection_times[dest] = row['EndTime']

    if infection_edges is not None:
        for n, e in infection_source_edge.items():
            infection_edges.append(e)

    # extract values of the dict
    timeline_number_of_infected = list(infection_times.values())
    timeline_number_of_infected.sort()

    return timeline_number_of_infected, infection_times


def compute_avg_prevalence(timeline, step, n_step, initial_time, n_nodes):
    """
    Computes the average prevalence for each iteration, i.e. simulation of infection

    :param timeline: list of results of timelines for the iterations
    :param step: lenght of each step
    :param n_step: number of steps for plotting
    :param initial_time: time of first flight's departure
    :param n_nodes: number of nodes in the network = number of airports

    :return: avg_timelines: list of n_step with number of cities infected for each step
    """

    sum_timelines = [0] * n_step

    for time in timeline:
        i_th = int(((time - initial_time)/step))
        sum_timelines[i_th] += 1

    avg_timelines = [x/n_nodes for x in np.cumsum(sum_timelines)]

    return avg_timelines


def simulate_multiple_iterations(n_iter, probs, seed_ids, first_dep, last_arr, n_nodes, n_step,
                                 infection_duration=-1, model='SI'):

    """
    Computes the vectors of data obtained by multiple simulations with different probabilities and seeds

    :param n_iter: number of iterations to be performed
    :param probs: one or multiple probabilities to be tested
    :param seed_ids: one or multiple seed_ids to be tested
    :param first_dep: time of first flight's departure
    :param last_arr: time of last flight's arrival
    :param n_nodes: number of nodes in the network
    :param n_step: number of steps in linspace

    for BONUS TASK
    :param model: type of model implemented
    :param infection_duration: duration of the infection - number of seconds

    :return datavecs: vectors of data to be plotted
    :return max_x: max x-value to be plotted, given by the difference between last_arr and firts_dep
    """

    step = int((last_arr-first_dep)/n_step)
    datavecs = []

    for p in probs:
        for seed_id in seed_ids:
            timelines = []
            for i in range(n_iter):
                if infection_duration == -1 and model == 'SI':
                    timeline, infection_times = simulate_infection(seed_id, p)
                else:
                    timeline, infection_times = simulate_infection(seed_id, p, infection_duration=infection_duration,
                                                                   model=model)
                avg_prevalence = compute_avg_prevalence(timeline, step, n_step, first_dep, n_nodes)
                timelines.append(avg_prevalence)
            # compute the average of each value over the iterations
            datavecs.append(list(np.mean(timelines, axis=0)))

    diff = last_arr-first_dep

    return datavecs, diff


def simulate_multiple_iterations_immunized(p, seed_ids, first_dep, last_arr, n_nodes, n_step, immu_set):
    """

    :param p: probability of infection
    :param seed_ids: list of int
    :param first_dep: time of first departure
    :param last_arr: time of last arrival
    :param n_nodes: number of nodes
    :param n_step: number of steps
    :param immu_set: set of immunized nodes
    :return: datavec: list of avg prevalence over the various iterations
    """
    step = int((last_arr-first_dep)/n_step)
    datavec, timelines = [], []

    for seed_id in seed_ids:
        timeline, infection_times = simulate_infection(seed_id, p, immunized=immu_set)
        avg_prevalence = compute_avg_prevalence(timeline, step, n_step, first_dep, n_nodes)
        timelines.append(avg_prevalence)

    datavec = list(np.mean(timelines, axis=0))
    return datavec


def get_network_measures(net, immunized=None):
    """
    Computes measures needed for task 4

    :param net: network
    :param immunized: int of number of nodes to be immunized
    :return: c: list of unweighted clustering coefficient for each node
             k: list of degrees for each node
             s: list of strengths for each node
             bc: list of betweenness centralities for each node

    """

    c, k, s, bc = [], [], [], []
    clustering_coeff = nx.clustering(net)
    # print(clustering_coeff)
    degrees = nx.degree(net)
    strengths = nx.degree(net, weight='weight')
    betweenness_centr = nx.betweenness_centrality(net, normalized=True)
    if immunized is None:
        for n in nodes:
            c.append(clustering_coeff[n])
            k.append(degrees[n])
            s.append(strengths[n])
            bc.append(betweenness_centr[n])

    else:
        ordered_c = sorted(clustering_coeff.items(), key=lambda x: x[1], reverse=True)
        c = set([int(x[0]) for x in ordered_c[:immunized]])
        ordered_k = sorted(degrees, key=lambda x: x[1], reverse=True)
        k = set([int(x[0]) for x in ordered_k[:immunized]])
        ordered_s = sorted(strengths, key=lambda x: x[1], reverse=True)
        s = set([int(x[0]) for x in ordered_s[:immunized]])
        ordered_bc = sorted(betweenness_centr.items(), key=lambda x: x[1], reverse=True)
        bc = set([int(x[0]) for x in ordered_bc[:immunized]])

    return c, k, s, bc


def get_link_measures(net):
    """
    Compute weights and edges betweenness centralities

    :param net:  network
    :return: w: list of weights
             eb: list of edge betweeenness centralities
    """
    w, eb, eb_w, eb_w2, eb_pr, eb_cl, eb_ev, eb_s = [], [], [], [], [], [], [], []

    # Edge weight and unweighted betweenness centrality
    edges = net.edges(data=True)
    betweenness_centr = nx.edge_betweenness_centrality(net, normalized=True)
    for e in edges:
        w.append(e[2]['weight'])
        eb.append(betweenness_centr[(e[0], e[1])])

    # Create a copy of the graph with inverse weights, square root is used to reduce the impact of high weights
    net1 = net.copy()
    edges1 = net1.edges(data=True)

    for e in edges1:
        w_e = e[2]['weight']
        net1[e[0]][e[1]]['weight'] = 1/(w_e**(1/3))

    # Weighted betweenness centrality on net1
    betweenness_centr_w = nx.edge_betweenness_centrality(net1, normalized=True, weight='weight')

    for e in edges1:
        eb_w.append(betweenness_centr_w[(e[0], e[1])])

    # Node dictionary for k-shells
    dict_k_shell = {}
    max_degree = max([net.degree(n) for n in net.nodes])

    for k in reversed(range(max_degree+1)):
        k_shell = nx.k_shell(net1, k=k)
        k_shell_nodes = k_shell.nodes()
        for i in k_shell_nodes:
            if i not in dict_k_shell:
                dict_k_shell[i] = k

    # node dict for pagerank
    dict_page_rank = nx.pagerank(net1, weight='weight')

    # closeness centrality
    closeness_centr = nx.closeness_centrality(net,  distance='weight')
    closeness_centr = dict(sorted(closeness_centr.items(), key=lambda pair: list(nodes).index(pair[0])))
    # eigenvector centrality
    eigenvector_centr = nx.eigenvector_centrality(net, tol=10**-1, weight='weight')
    eigenvector_centr = dict(sorted(eigenvector_centr.items(), key=lambda pair: list(nodes).index(pair[0])))
    # strengths of nodes
    strengths = dict(nx.degree(net1, weight='weight'))

    # For each edge, take lower value of centrality measureof the two nodes and use it to normalize previously
    # computed weighted betweenness
    j = 0
    for e in edges:
        eb_w2.append(eb_w[j]/min(dict_k_shell[e[0]], dict_k_shell[e[1]]))
        eb_pr.append(eb_w[j]/min(dict_page_rank[e[0]], dict_page_rank[e[1]]))
        eb_cl.append(eb_w[j]/min(closeness_centr[e[0]], closeness_centr[e[1]]))
        eb_ev.append(eb_w[j]/min(eigenvector_centr[e[0]], eigenvector_centr[e[1]]))
        eb_s.append(eb_w[j]/min(strengths[e[0]], strengths[e[1]]))
        j = j+1

    return w, eb, eb_w, eb_w2, eb_pr, eb_cl, eb_ev, eb_s


def compute_infection_medians(seed_ids, p, nodes):
    """
    Compute medians for each node over several simulation computed with random seed ids

    :param seed_ids: list of random seed ids
    :param p: probability of a node to be infected
    :param nodes: airports/nodes of the network
    :return: infections: list of medians
    """

    infections = {}
    for n in nodes:
        infections[n] = []

    for seed_id in seed_ids:
        timeline_number_of_infected, infection_times = simulate_infection(seed_id, p)
        for k, v in infection_times.items():
            infections[str(k)].append(v)

    res = []
    for n in infections:
        res.append(np.median(infections[n]))

    # list(infections.values())
    return res


def clean_values(infection_medians, c, k, s, bc):

    """
    Cleans the arrays from nodes values that have not been infected.
    :param infection_medians: list of median times of infections
    :param c: list of unweighted clustering coefficient for each node
    :param k: list of degrees for each node
    :param s: list of strengths for each node
    :param bc: list of betweenness centralities for each node
    :return: x_values: list of numpy arrays
             y_values: numpy array of median times of infections
    """

    y, c_new, k_new, s_new, bc_new = [], [], [], [], []
    cnt = 0
    for value in infection_medians:
        if not np.isnan(value):
            y.append(value-first_dep)
            c_new.append(c[cnt])
            k_new.append(k[cnt])
            s_new.append(s[cnt])
            bc_new.append(bc[cnt])
        cnt += 1

    y_values = np.array(y)
    x_values = [np.array(c_new), np.array(k_new), np.array(s_new), np.array(bc_new)]

    return x_values, y_values


def get_immunized_nodes(strategies, n_immu, net):
    """
    Computes the set

    :param strategies: list of strings - types of immunization strategy
    :param n_immu: number of nodes that are in the set of immunized nodes
    :param net: network

    :return: immunized_sets: list of sets of int representing node_ids
             immunized_nodes: set of int - total immunized nodes
    """
    immunized_sets, immunized_nodes = [], set()
    social, random, c, k, s, bc = set(), set(), set(), set(), set(), set()

    for strategy in strategies:
        if strategy == 'social_network':
            for i in range(n_immu):
                # take a random node from the network
                n = rd.choice(list(net.nodes()))
                neighbors = list(net.neighbors(n))
                # pick a random neighbor to be immunized
                immune = int(rd.choice(neighbors))
                social.add(immune)
        elif strategy == 'random':
            for i in range(n_immu):
                # take a random node form the network
                immune = int(rd.choice(list(net.nodes())))
                random.add(immune)
        else:
            c, k, s, bc = get_network_measures(net, n_immu)

    list_sets = [social, random, c, k, s, bc]
    for i_set in list_sets:
        immunized_sets.append(i_set)
        immunized_nodes = immunized_nodes.union(i_set)

    return immunized_sets, immunized_nodes


# ========================== MAIN ===========================

if __name__ == '__main__':

    network_path = 'aggregated_US_air_traffic_network_undir.edg'
    events_file = 'events_US_air_traffic_GMT.txt'
    airport_info = 'US_airport_id_info.csv'
    bg_figname = 'US_air_bg.png'

    net = nx.read_weighted_edgelist(network_path)
    # noinspection PyTypeChecker
    id_data = np.genfromtxt(airport_info, delimiter=',', dtype=None, names=True, encoding=None)
    xycoords = {}
    for row in id_data:
        xycoords[str(row['id'])] = (row['xcoordviz'], row['ycoordviz'])
    net = nx.read_weighted_edgelist(network_path)

    # noinspection PyTypeChecker
    event_data = np.genfromtxt(events_file, names=True, dtype=int)
    # sort all flights in increasing order of departure time
    event_data.sort(order=['StartTime'])

    # some parameters needed after multiple times
    first_dep = event_data[0]['StartTime']
    print("The first departure time is "+str(first_dep))
    last_arr = event_data[len(event_data)-1]['EndTime']
    print("The last arrival time is "+str(last_arr))
    print("The delta is "+str(last_arr-first_dep))
    nodes = list(net.nodes())
    edges = list(net.edges())
    n_nodes = len(nodes)

    ''' Task 1 '''
    # a) Find out at which time node_id=4 is infected if seed_id=0
    timeline_number_of_infected, infection_times = simulate_infection(seed_id=0, p=1.0)
    print('Task 1 (a) Anchorage is infected at time '+str(infection_times[41]))

    ''' Task 2 '''
    n_step = 100
    probs = [0.01, 0.05, 0.1, 0.5, 1.0]
    datavecs, max_x = simulate_multiple_iterations(n_iter=10, probs=probs, seed_ids=[0], first_dep=first_dep,
                                                   last_arr=last_arr, n_nodes=n_nodes, n_step=n_step)

    # a) plot fraction of infected nodes as a function of time for each of the infection probabilities 5 curves in one graph

    name = 'Averaged_prevalence'
    labels = ['p = 0.01', 'p = 0.05', 'p = 0.1', 'p = 0.5', 'p = 1.0']
    xlabel = 'time'
    ylabel = 'averaged prevalence'
    path = './averaged_prevalence_prob.png'
    colors = ['r', 'b', 'g', 'c', 'm']

    print("Task 2 (a) Plotting averaged prevalence over different probabilities with seed node Allentown (id=0)")

    fig1 = plot_averaged_prevalence(datavecs, labels, xlabel, ylabel, name, max_x, n_step, colors)
    fig1.savefig(path)
    plt.close(fig1)

    ''' Task 3 '''
    # a)
    seed_ids = [0, 4, 41, 100, 200]
    datavecs, max_x = simulate_multiple_iterations(n_iter=10, probs=[0.1], seed_ids=seed_ids, first_dep=first_dep,
                                                   last_arr=last_arr, n_nodes=n_nodes, n_step=n_step)

    name = 'Averaged_prevalence'
    labels = ['seed_id = 0', 'seed_id = 4', 'seed_id = 41', 'seed_id = 100', 'seed_id = 200']
    xlabel = 'time'
    ylabel = 'averaged prevalence'
    path = './averaged_prevalence_seed.png'
    colors = ['r', 'b', 'g', 'c', 'm']

    print("Task 3 (a) Plotting  averaged prevalence over different seed nodes with probability of infection 0.1")

    fig2 = plot_averaged_prevalence(datavecs, labels, xlabel, ylabel, name, max_x, n_step, colors)
    fig2.savefig(path)
    plt.close(fig2)

    ''' Task 4 '''
    # selecting a refuge
    p = 0.5
    n_iter = 50
    seed_ids = [int(rd.choice(nodes)) for i in range(n_iter)]

    # a)
    infection_medians = compute_infection_medians(seed_ids, p, nodes)
    c, k, s, bc = get_network_measures(net)

    x_values, y_values = clean_values(infection_medians, c, k, s, bc)

    prefix = 'Scatter of median times vs '
    names = [prefix+'clustering coefficient', prefix+'degree', prefix+'strength', prefix+'betweenness centrality']
    xlabels = ['unweighted clustering coefficient', 'degree', 'strength', 'unweighted betweenness centrality']
    ylabel = 'median infection time'
    pref = './scatter_median_'
    paths = [pref+'cluster_coeff.png', pref+'degree.png', pref+'strength.png', pref+'betw_centr.png']
    colors = ['r', 'g', 'b', 'c']

    figures = create_scatter_plots(names=names, x_values=x_values, y_values=y_values, colors=colors,
                                   x_labels=xlabels, y_label=ylabel)
    for fig, path in zip(figures, paths):
        fig.savefig(path)
        # plt.close(fig)

    # b) best predictor for infection times
    max_corr, best_pred = 0, ""
    for x_val, type_pred in zip(x_values, xlabels):
        corr, p_val = spearmanr(x_val, y_values)
        print(corr)
        if abs(corr) > abs(max_corr):
            max_corr = corr
            best_pred = type_pred

    print("Task 4 (b) The best predictor for infection times is " + best_pred + " with correlation " + str(max_corr))

    ''' Task 5 '''
    # Shutting down airports
    # immunization strategies:
    # 1. social network = pick a random node from the network and immunize a random neighbour of this node
    # 2. random node = immunize a random node
    # 3. unweighted clustering coefficient = pick nodes with higher value
    # 4. degree = pick nodes with higher value
    # 5. strength = pick nodes with higher value
    # 6. unweighted betweenness centrality = pick nodes with higher value

    # a) For each immunization strategy simulate the infection n_iter times and plot the averaged prevalence as a
    # function of time
    p = 0.5
    n_seeds = 20
    n_immu = 10
    seed_ids = []

    immunization_strategies = ['social_network', 'random', 'clustering_coeff', 'degree', 'strength', 'betweenness_centr']

    # compute the set of immune airports
    immunized_sets, immunized_nodes = get_immunized_nodes(immunization_strategies, n_immu, net)
    # calculate the 20 seed nodes, so that they dont belong to any of the immunized ones
    for i in range(n_seeds):
        n = int(rd.choice(nodes))
        while n in immunized_nodes:
            n = int(rd.choice(nodes))
        seed_ids.append(n)

    datavecs = []
    for immu_set in immunized_sets:
        # simulate infection multiple times
        datavec = simulate_multiple_iterations_immunized(p=p, seed_ids=seed_ids, first_dep=first_dep, last_arr=last_arr,
                                                         n_nodes=n_nodes, n_step=n_step, immu_set=immu_set)
        datavecs.append(datavec)

    max_x = last_arr-first_dep
    # plot averaged prevalences of the disease as a function of time
    name = 'Averaged_prevalence'
    labels = immunization_strategies
    xlabel = 'time'
    ylabel = 'averaged prevalence'
    path = './averaged_prevalence_immunization_strategies.png'
    colors = ['r', 'b', 'g', 'c', 'm', 'y']

    print("Task 5 (a) Plotting  averaged prevalence over different immunization strategies with probability of infection 0.5")

    fig4 = plot_averaged_prevalence(datavecs, labels, xlabel, ylabel, name, max_x, n_step, colors)
    fig4.savefig(path)
    plt.close(fig4)

    ''' Task 6 '''
    n_seeds = 20
    p = 0.5
    seed_ids = [int(rd.choice(nodes)) for i in range(n_seeds)]

    infection_edges = {}  # dict with key = (airport_a, airport_b) and value = # of iterations infection passes through
    # link out of all simulations performed
    # always a < b so that each link is considered one time for both directions

    # For each simulation, record which links are used to infect yet uninfected airports
    # (either susceptible airports or the infecting flight arrives before the already recorded infection time).
    for seed_id in seed_ids:
        infection_edges_single_iteration = []
        simulate_infection(seed_id=seed_id, p=p, infection_edges=infection_edges_single_iteration)
        for k in infection_edges_single_iteration:
            if k in infection_edges:
                infection_edges[k] += 1
            else:
                infection_edges[k] = 1

    # a) compute the fraction of times that each link is used for infecting the disease f_ij
    # first compute the total number of infection

    fractions = {}
    for k, v in infection_edges.items():
        key = (str(k[0]), str(k[1]))
        fractions[key] = v/n_seeds

    linewidths = []
    for e in edges:
        if e in fractions:
            linewidths.append(fractions[e])
        else:
            linewidths.append(0)

    # plot the network
    fig, ax = plot_network_usa(net, xycoords, bg_figname=bg_figname, edges=edges, linewidths=linewidths)
    fig.savefig('./US_air_network.png')

    # c) Create scatter plots showing f_ij as a function of the following link properties:
    #       i) link weight w_ij
    #       ii) unweighted link betweenness centrality eb_ij
    w, eb, eb_w, eb_w2, eb_pr, eb_cl, eb_ev, eb_s = get_link_measures(net)

    x_values = [np.array(w), np.array(eb)]
    y_values = np.array(linewidths)

    prefix = 'Scatter of fraction of infecting times vs '
    names = [prefix+'link weight', prefix+'link betweenness centrality']
    xlabels = ['weight', 'unweighted link betweenness centrality']
    ylabel = 'fraction of infecting times'
    pref = './scatter_fraction_'
    paths = [pref+'link_weight.png', pref+'link_betw_centr.png']
    colors = ['r', 'g']

    figures = create_scatter_plots(names=names, x_values=x_values, y_values=y_values, colors=colors,
                                   x_labels=xlabels, y_label=ylabel)
    for fig, path in zip(figures, paths):
        fig.savefig(path)
        # plt.close(fig)

    # Compute the spearman correlation coefficient between f_ij and the two link-wise measures
    max_corr, best_pred = 0, ""
    for x_val, type_pred in zip(x_values, xlabels):
        corr, p_val = spearmanr(x_val, y_values)
        print(corr)
        if abs(corr) > abs(max_corr):
            max_corr = corr
            best_pred = type_pred

    print("Task 6 (c) The best predictor for fraction f_ij is " + best_pred + " with correlation " + str(max_corr))

    ''' BONUS Task 1'''

    x_values = [np.array(eb_w), np.array(eb_w2), np.array(eb_pr), np.array(eb_cl), np.array(eb_ev), np.array(eb_s)]
    y_values = np.array(linewidths)

    xlabels = ['weighted link betweenness centrality', 'k-shell', 'page_rank', 'closeness', 'eigenvector', 'strength']
    names = [prefix+'weighted_bc', prefix+'k-shell', prefix+'page_rank', prefix+'closeness', prefix+'eigenvector',
             prefix+'strenth']
    paths = [pref+'weighted_bc.png', pref+'k-shell.png', pref+'page_rank.png', pref+'closeness.png', pref+'eigenvector.png',
             pref+'strenth.png']
    colors = ['r', 'g', 'b', 'c', 'm', 'y']

    figures = create_scatter_plots(names=names, x_values=x_values, y_values=y_values, colors=colors,
                                   x_labels=xlabels, y_label=ylabel)
    for fig, path in zip(figures, paths):
        fig.savefig(path)

    # Compute the spearman correlation coefficient between f_ij and the two link-wise measures
    max_corr, best_pred = 0, ""
    for x_val, type_pred in zip(x_values, xlabels):
        corr, p_val = spearmanr(x_val, y_values)
        print(corr)
        if abs(corr) > abs(max_corr):
            max_corr = corr
            best_pred = type_pred

    print("BONUS TASK 1 - The best predictor for fraction f_ij is " + best_pred + " with correlation " + str(max_corr))

    ''' BONUS TASK 2 '''

    infection_durations = [15000, 50000, 100000, 200000]  # 4h ish, 1/2 day 13h ish, more than a day 27h, more than 2 days 55h
    models = ['SIS', 'SIR']

    for model in models:
        for infection_duration in infection_durations:
            seed_ids = [0, 4, 41, 100, 200]
            datavecs, max_x = simulate_multiple_iterations(n_iter=10, probs=[0.5], seed_ids=seed_ids, first_dep=first_dep,
                                                           last_arr=last_arr, n_nodes=n_nodes, n_step=n_step,
                                                           infection_duration=infection_duration, model=model)

            name = 'Averaged_prevalence'
            labels = ['seed_id = 0', 'seed_id = 4', 'seed_id = 41', 'seed_id = 100', 'seed_id = 200']
            xlabel = 'time'
            ylabel = 'averaged prevalence'
            path = './BONUS_averaged_prevalence_seed_'+model+'_'+str(infection_duration)+'.png'
            colors = ['r', 'b', 'g', 'c', 'm']

            print("BONUS TASK Plotting  averaged prevalence for " + model + " model with infection duration of " +
                  str(infection_duration))

            fig = plot_averaged_prevalence(datavecs, labels, xlabel, ylabel, name, max_x, n_step, colors)
            fig.savefig(path)
            plt.close(fig)

import networkx as nx
import matplotlib.pyplot as plt
import math
import sys


n_tops = 10000
min_pval = 0.001
min_enrich = 2
min_signif = 5
filter_parent = True
min_degree = 2
topGO_file = 'mol_fun/topGO_results_fisher'
filter_file = ''
printit = True
size_corr = 5
font_size = 8
iterations = 10000
n_plots = 1
empty_nodes = 0

class go_explore:

    def __init__(self, topGO, connections):
        self.topGO = topGO
        self.connections = connections
        self.indentation = ''
        
    def get_signif_parents(self, go):
        top_data = self.topGO[go]
        to_return = []
        if go not in self.connections:
            return
        for parent in self.connections[go]:
            if parent not in self.topGO:
                continue
            parent_data = self.topGO[parent]
            if filter_parent and (parent_data['pvalue'] > min_pval or parent_data['enrichment'] < min_enrich or parent_data['significant'] < min_signif):
                continue
            to_return.append([go, self.connections[go][parent], parent])
            for x in self.get_signif_parents(parent):
                to_return.append(x)
        return to_return 


def match_split(line, start):
    if linea.startswith(start):
        return linea.lstrip(start).split(' ')[0]

go_connections = {}
info = {'id: ':'', 'is_a: ':'is_a', 'relationship: part_of ':'is_part_of', 'relationship: regulates ':'regulates', 'relationship: negatively_regulates ':'regulates', 'relationship: positively_regulates ':'regulates'}
go = ''
for linea in open('/Users/martin/notDropbox/utils/goTerms/go.obo'):
    linea = linea.rstrip('\n')
    for k in info:
        match_info = match_split(linea, k)
        if not match_info:
            continue
        if k != 'id: ':
            if go not in go_connections:
                go_connections[go] = {}
            go_connections[go][match_info] = info[k]
        else:
            go = match_info
        break

n = 0
tops = []
topGO_data = {}
fields = ['significant', 'enrichment', 'pvalue', 'name']
for linea in open(topGO_file):
    if linea.startswith('GO.ID'):
        continue
    linea = linea.rstrip('\n')
    splat = linea.split('\t')
    topGO_data[splat[0]] = dict(zip(fields,[float(splat[3]), float(splat[3])/float(splat[4]), float(splat[5]), splat[1]])) #.replace(' ', '_')]))
    if n > n_tops or float(splat[5]) > min_pval or float(splat[3])/float(splat[4]) < min_enrich or float(splat[3]) < min_signif:
        continue
    tops.append(splat[0])
    n += 1

filter_out = []
if filter_file:
    for linea in open(filter_file):
        linea = linea.rstrip('\n')
        filter_out.append(linea.split(' ')[0])

G = nx.Graph()
go_exp = go_explore(topGO_data, go_connections)
for k in tops:
    connections = go_exp.get_signif_parents(k)
    if not connections:
        continue
    for (node1, connection, node2) in connections:
        G.add_edge(node1, node2)
        G[node1][node2]['connection'] = connection

G2 = nx.Graph()
degree = G.degree()
for edge in G.edges():
    node1, node2 = edge
    if degree[node1] >= min_degree or degree[node2] >= min_degree:
        if node1 in filter_out or node2 in filter_out:
            continue
        G2.add_edge(node1, node2)

for n in range(empty_nodes):
    G2.add_node(n)
    topGO_data[n] = {'name':'', 'significant':0, 'pvalue':0.5}

node_colors, sizes, labels = [], [], {}
for node in G2.nodes():
    labels[node] = topGO_data[node]['name']
    sizes.append(topGO_data[node]['significant'])
    node_colors.append(-math.log10(topGO_data[node]['pvalue']))

edges = G2.edges()
edge_colors = ['black'] * len(edges)
for n in range(len(edges)):
    if G[edges[n][0]][edges[n][1]]['connection'] == 'regulates':
        edge_colors[n] = 'red'

if printit:
    for label in labels:
        print label, labels[label] 

for k in range(n_plots):
    nx.draw(G2, edge_color=edge_colors, labels=labels, node_size=[x*size_corr for x in sizes], node_color=node_colors, cmap=plt.cm.Blues, font_size=font_size, iterations=iterations)
    plt.savefig('path' + str(k) + '.svg')
    plt.close()

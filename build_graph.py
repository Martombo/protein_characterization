n_tops = 10000
min_pval = 0.05
min_enrich = 2
min_signif = 3
filter_parent = True
lone_node = False

class go_explore:

    def __init__(self, topGO, connections):
        self.topGO = topGO
        self.connections = connections
        self.indentation = ''
        
    def get_signif_parents(self, go):
        top_data = self.topGO[go]
        to_print = self.indentation + top_data['name'] + ' ' + str(top_data['pvalue']) + ' ' + str(top_data['significant']) + '\n'
        self.indentation += '--'
        if go not in self.connections:
            return
        for parent in self.connections[go]:
            if parent not in self.topGO:
                continue
            parent_data = self.topGO[parent]
            if filter_parent and (parent_data['pvalue'] > min_pval or parent_data['enrichment'] < min_enrich or parent_data['significant'] < min_signif):
                continue
            to_print += self.get_signif_parents(parent)
        self.indentation = self.indentation[2:]
        return to_print 

    def get_signif_parents_sif(self, go):
        top_data = self.topGO[go]
        to_print = ''
        if lone_node:
            to_print = top_data['name'] + '\n'
        if go not in self.connections:
            return
        for parent in self.connections[go]:
            if parent not in self.topGO:
                continue
            parent_data = self.topGO[parent]
            if filter_parent and (parent_data['pvalue'] > min_pval or parent_data['enrichment'] < min_enrich or parent_data['significant'] < min_signif):
                continue
            to_print += top_data['name'] + ' ' + self.connections[go][parent] + ' ' + parent_data['name'] + '\n'
            to_print += self.get_signif_parents_sif(parent)
        return to_print 


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
for linea in open('top_w01F'):
    if linea.startswith('GO.ID'):
        continue
    linea = linea.rstrip('\n')
    splat = linea.split('\t')
    topGO_data[splat[0]] = dict(zip(fields,[int(splat[3]), int(splat[3])/float(splat[4]), float(splat[5]), splat[1].replace(' ', '_')]))
    if n > n_tops or float(splat[5]) > min_pval or int(splat[3])/float(splat[4]) < min_enrich or int(splat[3]) < min_signif:
        continue
    tops.append(splat[0])
    n += 1

go_exp = go_explore(topGO_data, go_connections)
for k in tops:
    go_parents = go_exp.get_signif_parents_sif(k)
    if go_parents:
        print(go_parents)

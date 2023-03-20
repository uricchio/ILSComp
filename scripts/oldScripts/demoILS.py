import dendropy
import random
import numpy as np
from scipy.stats import pearsonr
from dendropy.simulate import treesim
from dendropy.calculate import treecompare
from dendropy.datamodel import taxonmodel 
from dendropy.model import reconcile

# times are in units of N generations (i.e., 2 is 2N generations, or one coalescent unit)
#sp_tree_str = """\
#[&R] (A:10,(B:6,(C:4,(D:2,(E:1,F:1):2):2):2):4);
#"""

sp_tree_str = """\
[&R] ('Penicillium-cvjetkovicii-v1':0.223847,(((('Penicillium-commune-v1':0.03886,'Penicillium-solitum-v1':0.032132):0.005494,(('Penicillium-mb-v1':0.013052,'Penicillium-polonicum-v1':0.005086):0.021175,'Penicillium-verrucosum-v1':0.040045):0.004422):0.033079,'Penicillium-chrysogenum-v1':0.079855):0.1291,'Penicillium-brevicompactum-v1':0.266037):0.223847); 
"""

# set population size, which we need to get the units right
sp_tree = dendropy.Tree.get(data=sp_tree_str, schema="newick")
for edge in sp_tree.postorder_edge_iter():
        edge.pop_size = 0.5

# make the gene tree map
mapDict = {}
for node in sp_tree.leaf_node_iter():
    mapDict[node.taxon] = node.taxon

myMap = taxonmodel.TaxonNamespaceMapping(mapping_dict=mapDict)

# make the gene tree
gene_tree = treesim.contained_coalescent_tree(containing_tree=sp_tree,gene_to_containing_taxon_map=myMap)
sp_tree.migrate_taxon_namespace(gene_tree.taxon_namespace)

rf = treecompare.symmetric_difference(gene_tree, sp_tree)

# if you want to print out the tree, uncomment below
#print(gene_tree.as_ascii_plot())

# function to map num of seq diffs to tree
# adapted from https://dendropy.org/primer/trees.html#the-tree-class
# basically we are simulating a quant trait that drifts on tree
def evolve_char(tree, theta=200, mu = 0, sd = 1):
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node.value = 0
        else:
            numVar = np.random.poisson(theta*node.edge.length)
            node.value = node.parent_node.value+random.gauss(numVar*mu,(numVar**0.5)*sd)
    return tree

# run the sim on the tree
gene_tree = evolve_char(gene_tree)
pdmGT = gene_tree.phylogenetic_distance_matrix()
pdmST = sp_tree.phylogenetic_distance_matrix()

STdists = []
GTdists = []
trDist = []

for node in gene_tree.leaf_node_iter():
    for node2 in gene_tree.leaf_node_iter():
        if node2.taxon > node.taxon: # or node2.taxon == node.taxon:
            STdists.append( pdmST.distance(myMap[node.taxon], myMap[node2.taxon]))
            GTdists.append( pdmGT.distance(myMap[node.taxon], myMap[node2.taxon]))
            trDist.append(1+abs(node.value-node2.value))
            #print(myMap[node.taxon], myMap[node2.taxon], pdmST.distance(myMap[node.taxon], myMap[node2.taxon]),pdmGT.distance(node.taxon,node2.taxon), 1+abs(node.value-node2.value), gene_tree.as_string(schema="newick")[:-1], rf)

for i in range(len(STdists)):
    print (STdists[i], GTdists[i], trDist[i])

print ("#", pearsonr(STdists, trDist), pearsonr(GTdists, trDist), rf)
print( "#GENE: ", gene_tree.as_string(schema="newick")[:-1], rf)
print( "#SPECIES: ", sp_tree.as_string(schema="newick")[:-1], rf)


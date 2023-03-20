import dendropy
import random
import sys
import numpy as np

from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import poisson
from dendropy.simulate import treesim
from dendropy.calculate import treecompare
from dendropy.datamodel import taxonmodel 
from dendropy.model import reconcile
from dendropy.simulate import treesim

# helper function
def comp(bip):
    ret = ''
    for c in bip:
        if c == '0':
            ret += '1'
        if c == '1':
            ret += '0'
    return ret

#inputs
numSims = int(sys.argv[1])
pop_size = float(sys.argv[2])
theta = 0.1/pop_size

# times are in units of N generations (i.e., 2 is 2N generations, or one coalescent unit)
sp_tree_str = """\
[&R] (Penicillium-brevicompactum:0.1497010058,(Penicillium-chrysogenum:0.0516743022,((Penicillium-commune:0.0132466035,Penicillium-solitum:0.0149358607):0.0041523096,((Penicillium-mb:0.0093122349,Penicillium-polonicum:0.0074256451):0.0103630038,Penicillium-verrucosum:0.0295819813):0.0020373155):0.0275815073):0.0745440008,Penicillium-cvjetkovicii:0.3495454216);
"""

# set population size, which we need to get the units right
sp_tree = dendropy.Tree.get(data=sp_tree_str, schema="newick")
for edge in sp_tree.postorder_edge_iter():
    edge.pop_size = pop_size

# uncomment to simulate a random tree
#myNames = ["z"+str(i) for i in range(8)]
#taxa = dendropy.TaxonNamespace(myNames)
#tree = treesim.pure_kingman_tree(
#        taxon_namespace=taxa,
#        pop_size=1)
#sp_tree = treesim.birth_death_tree(birth_rate=0.1, death_rate=0.1, num_extant_tips=8)
#sp_tree.as_string(schema='newick').split()[1]

# make the gene tree map
mapDict = {}
for node in sp_tree.leaf_node_iter():
    mapDict[node.taxon] = node.taxon

#print(sp_tree.as_ascii_plot())
#print(sp_tree.as_string(schema='newick').split()[1])

myMap = taxonmodel.TaxonNamespaceMapping(mapping_dict=mapDict)

# make the gene tree
gene_tree = treesim.contained_coalescent_tree(containing_tree=sp_tree,gene_to_containing_taxon_map=myMap)
sp_tree.migrate_taxon_namespace(gene_tree.taxon_namespace)

gene_tree.encode_bipartitions()
sp_tree.encode_bipartitions()

# make a set of bipartitions and a map of species names
bips = {}
specMap = {}
for key in gene_tree.bipartition_edge_map:
    bips[str(key)] = 0
    bips[comp(str(key))] = 0
    if gene_tree.bipartition_edge_map[key].head_node.taxon is not None:
        i = 0
        for c in str(key):
            if str(key)[i] == "1":
                break
            i += 1
        specMap[gene_tree.bipartition_edge_map[key].head_node.taxon] = i

rf = treecompare.symmetric_difference(gene_tree, sp_tree)

# if you want to print out the tree, uncomment below
#print(gene_tree.as_ascii_plot())

# function to map num of seq diffs to tree
# adapted from https://dendropy.org/primer/trees.html#the-tree-class
# basically we are simulating a quant trait that drifts on tree

def evolve_binary_trait(tree, theta=10/pop_size):
 
    tot = 0
    for node in tree.preorder_node_iter():
        tot += node.edge.length
     
    # get an appropriate p for the tree height
    totExp = 1
    p = totExp/(theta*tot)

    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node.value = 0
            node.numVar = 0
        else:
            node.numVar = np.random.poisson(theta*node.edge.length)
            numSwitches = np.random.poisson(node.numVar*p)
            if numSwitches % 2 == 0:
                node.value = node.parent_node.value
            else:
                node.value = 1-node.parent_node.value
            #print(numSwitches,node.taxon,node.value, node.numVar, node)
    return tree

def normalize(tree):
    vals = []
    for node in tree.leaf_node_iter():
        vals.append(node.value)
    m = np.mean(vals)
    s = np.std(vals)
    for node in tree.leaf_node_iter():
        node.value = (node.value-m)/s
    return tree

# run the sim on the tree
gene_tree = evolve_binary_trait(gene_tree, theta=theta)

# get bip of gene tree traits
gene_tree_tr_bip = [0 for i in gene_tree.leaf_node_iter()]
tot = 0
i = 0
for node in gene_tree.leaf_node_iter():
    if node.value == 1:
        gene_tree_tr_bip[specMap[node.taxon]] = 1
        tot += 1
    i += 1

# only keep non-trivial bips
if tot == 0 or tot == i or tot ==1 or tot == i-1:
    exit()

gene_tree_tr_bip_str = ''
for spec in gene_tree_tr_bip:
    gene_tree_tr_bip_str  += str(spec)

nsubs = 0
in_tree = 0
in_Stree = 0

# check if bip in bips
if gene_tree_tr_bip_str in bips or comp(str(key)) in bips:
    in_tree = 1
    for key in gene_tree.bipartition_edge_map:
        if str(key) == gene_tree_tr_bip_str or comp(str(key)) == gene_tree_tr_bip_str: 
            #print( gene_tree.bipartition_edge_map[key].head_node)
            nsubs += gene_tree.bipartition_edge_map[key].head_node.numVar
            break

print("#GENE_TREEs nsubs, bip, in_g_tree, rf:", nsubs, gene_tree_tr_bip_str, in_tree, rf)

for k in range(numSims):
    
    # pois param, i.e. the length of the bip
    tot_pois = 0
  
    # make new tree and encode bips
    gene_tree_new = treesim.contained_coalescent_tree(containing_tree=sp_tree,gene_to_containing_taxon_map=myMap)
    gene_tree_new = evolve_binary_trait(gene_tree_new)
    gene_tree_new.encode_bipartitions()

    # get bipartitions map of new gene tree
    bips = {}
    for key in gene_tree_new.bipartition_edge_map:
        bips[str(key)] = 0
        if gene_tree_new.bipartition_edge_map[key].head_node.taxon is not None:
            i = 0
            for c in str(key):
                if str(key)[i] == "1":
                    break
                i += 1

    # check if bip in bips
    if gene_tree_tr_bip_str in bips or comp(gene_tree_tr_bip_str) in bips:
        for key in gene_tree_new.bipartition_edge_map:
            if str(key) == gene_tree_tr_bip_str or comp(str(key)) == gene_tree_tr_bip_str:
                tot_pois += poisson.rvs(gene_tree_new.bipartition_edge_map[key].head_node.edge.length*theta)
                break
    

    rf = treecompare.symmetric_difference(gene_tree_new, gene_tree)
    print(tot_pois, rf)

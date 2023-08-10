import dendropy
import random
import sys
import numpy as np

from scipy.stats import pearsonr
from scipy.stats import spearmanr
from dendropy.simulate import treesim
from dendropy.calculate import treecompare
from dendropy.datamodel import taxonmodel 
from dendropy.model import reconcile
from dendropy.simulate import treesim

#inputs
numSims = int(sys.argv[1])
pop_size = float(sys.argv[2])

# times are in units of N generations (i.e., 2 is 2N generations, or one coalescent unit)
sp_tree_str = """\
[&R] (Penicillium-brevicompactum:0.1497010058,(Penicillium-chrysogenum:0.0516743022,((Penicillium-commune:0.0132466035,Penicillium-solitum:0.0149358607):0.0041523096,((Penicillium-mb:0.0093122349,Penicillium-polonicum:0.0074256451):0.0103630038,Penicillium-verrucosum:0.0295819813):0.0020373155):0.0275815073):0.0745440008,Penicillium-cvjetkovicii:0.3495454216);
"""


# set population size, which we need to get the units right
sp_tree = dendropy.Tree.get(data=sp_tree_str, schema="newick")
for edge in sp_tree.postorder_edge_iter():
        edge.pop_size = pop_size

# uncomment to simulate a random tree
#myNames = ["z"+str(i) for i in range(50)]
#taxa = dendropy.TaxonNamespace(myNames)
#tree = treesim.pure_kingman_tree(
#        taxon_namespace=taxa,
#        pop_size=1)
#sp_tree = treesim.birth_death_tree(birth_rate=0.1, death_rate=0.1, num_extant_tips=50)

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

rf = treecompare.symmetric_difference(gene_tree, sp_tree)

# if you want to print out the tree, uncomment below
#print(gene_tree.as_ascii_plot())

# function to map num of seq diffs to tree
# adapted from https://dendropy.org/primer/trees.html#the-tree-class
# basically we are simulating a quant trait that drifts on tree
def evolve_char(tree, theta=10/pop_size, mu = 0, sd = 1):
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node.value = 0
        else:
            numVar = np.random.poisson(theta*node.edge.length)
            node.value = node.parent_node.value+random.gauss(numVar*mu,(numVar**0.5)*sd)
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
gene_tree = evolve_char(gene_tree)
gene_tree = normalize(gene_tree)

pdmGT = gene_tree.phylogenetic_distance_matrix()
pdmST = sp_tree.phylogenetic_distance_matrix()


pears = []
STdists = []
for node in gene_tree.leaf_node_iter():
    GTdists = []
    trDist = []
    spec = []
    STdists = []
    for node2 in gene_tree.leaf_node_iter():
        if node2.taxon != node.taxon: # or node2.taxon == node.taxon:
            STdists.append( pdmST.distance(myMap[node.taxon], myMap[node2.taxon]))
            GTdists.append( pdmGT.distance(myMap[node.taxon], myMap[node2.taxon]))
            trDist.append((node.value-node2.value)**2)
            spec.append([node.taxon,node2.taxon]) 
    pears.append( [pearsonr(GTdists, trDist)[0], pearsonr(GTdists, trDist)[1]])
# print header
print ("#SPECIES_TREE:", pearsonr(STdists, trDist))
print( "#GENE_TREEs p:",  end=' ')
for thing in pears:
    print (thing[1], end=' ')
print() 
print( "#GENE_TREEs r:",  end=' ')
for thing in pears:
    print (thing[0], end=' ')
print()

fh = open("../simData/spec.txt", "w")
treeStr = sp_tree.as_string(schema='newick').split()[1]
treeFinal = ''
for c in treeStr:
    if c != "'":
        treeFinal += c
   
fh.write(treeFinal+"\n")
fh.close()

#for PGLS uncomment
#for i in range(len(STdists)):
#        print(STdists[i], GTdists[i], trDist[i], spec[i][0],spec[i][1], 0, rf)

for k in range(numSims):
    gene_tree_new = treesim.contained_coalescent_tree(containing_tree=sp_tree,gene_to_containing_taxon_map=myMap)

    rf = treecompare.symmetric_difference(gene_tree_new, sp_tree)
 
    pdmGT = gene_tree_new.phylogenetic_distance_matrix()
    #pdmST = sp_tree.phylogenetic_distance_matrix()

    for node in gene_tree.leaf_node_iter():
        STdists = []
        GTdists = []
        trDist = []
        spec = []

        for node2 in gene_tree.leaf_node_iter():
            if node2.taxon != node.taxon:
                STdists.append( pdmST.distance(myMap[node.taxon], myMap[node2.taxon]))
                GTdists.append( pdmGT.distance(myMap[node.taxon], myMap[node2.taxon]))
                trDist.append((node.value-node2.value)**2)
                spec.append([node.taxon,node2.taxon])          
  
        res = pearsonr(GTdists, trDist)
        print(rf, res[0], res[1])
    #for PGLS uncomment
    #for i in range(len(STdists)):
    #    print(STdists[i], GTdists[i], trDist[i], spec[i][0],spec[i][1], k+1, rf, pearsonr(GTdists, trDist))

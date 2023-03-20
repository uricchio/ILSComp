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

gene_tree.encode_bipartitions()

# make a set of bipartitions and a map of species names
bips = {}
specMap = {}
for key in gene_tree.bipartition_edge_map:
    bips[str(key)] = 0
    if gene_tree.bipartition_edge_map[key].head_node.taxon is not None:
        i = 0
        for c in str(key):
            if str(key)[i] == 1:
                break
            i += 1
        specMap[gene_tree.bipartition_edge_map[key].head_node.taxon] = i

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

def evolve_binary_trait(tree, theta=1/pop_size, p = 1):
 
    tot = 0
    for node in tree.preorder_node_iter():
        tot += node.edge.length
     
    # get an appropriate p for the tree height
    totExp = 10
    p = totExp/(theta*tot)

    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node.value = 0
            node.numVar = 0
        else:
            node.numVar = np.random.poisson(theta*node.edge.length)
            numSwitches = np.random.poisson(node.numVar*p)
            #print(numSwitches)
            if numSwitches % 2 == 0:
                node.value = node.parent_node.value
            else:
                node.value = 1-node.parent_node.value
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
gene_tree = evolve_binary_trait(gene_tree)

# get bip of gene tree traits
gene_tree_tr_bip = [0 for i in gene_tree.leaf_node_iter()]
i = 0
tot = 0
for node in gene_tree.leaf_node_iter():
    if node.value == 1:
        gene_tree_tr_bip[i] = 1
        tot += 1
    i += 1
if tot == 0 or tot == i or tot ==1 or tot == i-1:
    exit()
gene_tree_tr_bip_str = ''
for spec in gene_tree_tr_bip:
    gene_tree_tr_bip_str  += str(spec)


nsubs = 0
in_tree = 0
# check if bip in bips

if gene_tree_tr_bip_str in bips:
    for key in gene_tree.bipartition_edge_map:
        #print(str(key))
        if str(key) == gene_tree_tr_bip_str:
            in_tree = 1
            nsubs += gene_tree.bipartition_edge_map[key].head_node.numVar
            break
        
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

print( "#GENE_TREEs nsubs, bip, in_g_tree:", nsubs, gene_tree_tr_bip_str, in_tree)

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
    # make new tree and encode bips
    gene_tree_new = treesim.contained_coalescent_tree(containing_tree=sp_tree,gene_to_containing_taxon_map=myMap)
    gene_tree_new= evolve_binary_trait(gene_tree_new)
    gene_tree_new.encode_bipartitions()

    # get bipartitions map of new gene tree
    bips = {}
    for key in gene_tree_new.bipartition_edge_map:
        bips[str(key)] = 0
        if gene_tree_new.bipartition_edge_map[key].head_node.taxon is not None:
            i = 0
            for c in str(key):
                if str(key)[i] == 1:
                    break
                i += 1

    nsubs = 0
    # check if bip in bips
    if gene_tree_tr_bip_str in bips:
        for key in gene_tree_new.bipartition_edge_map:
            if str(key) == gene_tree_tr_bip_str:
                nsubs += gene_tree_new.bipartition_edge_map[key].head_node.numVar
                break

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
        print(rf, res[0], res[1], nsubs)
    #for PGLS uncomment
    #for i in range(len(STdists)):
    #    print(STdists[i], GTdists[i], trDist[i], spec[i][0],spec[i][1], k+1, rf, pearsonr(GTdists, trDist))

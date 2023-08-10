import random
import sys
import os
import numpy as np
import dendropy

from scipy.optimize import minimize
from scipy.stats import pearsonr
from scipy.stats import linregress
from scipy.stats import spearmanr
from scipy.stats import percentileofscore
from scipy.interpolate import interp1d
from dendropy.simulate import treesim
from dendropy.calculate import treecompare
from dendropy.datamodel import taxonmodel
from dendropy.model import reconcile
from dendropy.simulate import treesim
from collections import defaultdict


class SimulateILS():

    def __init__(self,popSize=100,numSims=100,theta=0.1):

        self.popSize = popSize
        self.numSims = numSims
        self.theta = theta
        self.sp_tree_str = """\
                           [&R] (Penicillium-cvjetkovicii:0.347413,(Penicillium-brevicompactum:0.149364,(Penicillium-chrysogenum:0.052808,((Penicillium-solitum:0.015629,Penicillium-commune:0.013995):0.004497,((Penicillium-mb:0.009725,Penicillium-polonicum:0.008057):0.010835,Penicillium-verrucosum:0.030704):0.002275):0.030322):0.076323):0.010049);\
                           """
        self.sp_tree = dendropy.Tree.get(data=self.sp_tree_str, schema="newick")       
        self.gene_tree = self.sp_tree
        self.gene_tree_list = []
        self.mapDict = {}
        self.corr = []
        self.real_corr = []
        self.nullCorr = {}
        self.pdmST = []
        self.vals = {}
        self.myTrees = []
        self.trDist = defaultdict(dict) 
        self.copyNum = []
        #self.reRootOutgroup()

    def reRootOutgroup(self,og='Penicillium-cvjetkovicii'):
        myOG = self.sp_tree.find_node_with_taxon_label(og)
        self.sp_tree.reroot_at_midpoint(update_bipartitions=False) #to_outgroup_position(myOG, update_bipartitions=False)   
        print(self.sp_tree.as_string(schema="newick"))   

    def setPopSize(self,popSize=-1):
        if popSize < 0:
            popSize = self.popSize
        for edge in self.sp_tree.postorder_edge_iter():
            edge.pop_size = popSize

    def initRandomTree(self, n=8):
        myNames = ["z"+str(i) for i in range(n)]
        taxa = dendropy.TaxonNamespace(myNames)
        tree = treesim.pure_kingman_tree(
            taxon_namespace=taxa,
            pop_size=1)
        sp_tree = treesim.birth_death_tree(birth_rate=0.1, death_rate=0.1, num_extant_tips=n)
        self.sp_tree = sp_tree

    def makeTaxonMap(self):
        mapDict = {}
        for node in self.sp_tree.leaf_node_iter():
            mapDict[node.taxon] = node.taxon
        self.mapDict = taxonmodel.TaxonNamespaceMapping(mapping_dict=mapDict)

    def simGeneTree(self):
        self.gene_tree = treesim.contained_coalescent_tree(containing_tree=self.sp_tree,gene_to_containing_taxon_map=self.mapDict)
        self.sp_tree.migrate_taxon_namespace(self.gene_tree.taxon_namespace)

    def calcRF(self,t1, t2):
        return treecompare.symmetric_difference(t1, t2)
        
    def outputCladogram(self, t):
     
        for edge in t.preorder_edge_iter():
            edge.length = None
        print(t.as_string(schema="newick")) 


    def evolveTrait(self, theta = 1, mu = 0, sd = 1):
        for node in self.gene_tree.preorder_node_iter():
            if node.parent_node is None:
                node.value = 0
            else:
                numVar = np.random.poisson(theta*node.edge.length)
                node.value = node.parent_node.value+random.gauss(numVar*mu,(numVar**0.5)*sd)

    def evolveBinaryTrait(self, tree, theta=1):

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

    def normalize(self, tree):
        vals = []
        for node in tree.leaf_node_iter():
            vals.append(node.value)
        m = np.mean(vals)
        s = np.std(vals)

        #randVar = np.random.normal(
        for node in tree.leaf_node_iter():
            node.value = (node.value-m)/s
        return tree

    def simCausalGene(self,pgls=False,geneSet = False):
        if not geneSet:
            self.evolveTrait(theta=self.theta,mu=0,sd=1)
            self.gene_tree = self.normalize(self.gene_tree)

        pdmGT = self.gene_tree.phylogenetic_distance_matrix()
        self.pdmST = self.sp_tree.phylogenetic_distance_matrix()

        corr = []
        STdists = []
        for node in self.gene_tree.leaf_node_iter():
            GTdists = []
            trDist = []
            spec = []
            STdists = []
            for node2 in self.gene_tree.leaf_node_iter():
                if node2.taxon.label != node.taxon.label: # or node2.taxon == node.taxon:
                    STdists.append( self.pdmST.distance(self.mapDict[node.taxon], self.mapDict[node2.taxon]))
                    GTdists.append( pdmGT.distance(self.mapDict[node.taxon], self.mapDict[node2.taxon]))
                    trDist.append((node.value-node2.value)**2)
                    spec.append([str(node.taxon),str(node2.taxon)])
            if pgls:
                #print("here")
                for i in range(len(STdists)):
                    self.vals[(spec[i][0],spec[i][1])] = [trDist[i], GTdists[i]]
                #print(len(self.vals))
            out = linregress(GTdists, trDist)
            corr.append( [out.rvalue, out.pvalue,str(node.taxon)[1:-1]])
        self.corr = corr

    def makeOutputHeader(self):
        # print header
        print( "#SPECIES:",  end=' ')
        for thing in self.corr:
            print (thing[2], end=' ')
        print()
        print( "#GENE_TREEs p:",  end=' ')
        for thing in self.corr:
            print (thing[1], end=' ')
        print()
        print( "#GENE_TREEs r:",  end=' ')
        for thing in self.corr:
            print (thing[0], end=' ')
        print()

    def makeOutputHeaderSet(self):
        print( "#GENE_TREEs p:",  end=' ')
        for thing in self.corr:
            print(thing[1], end= ' ')
        print()
        print( "#GENE_TREEs r:",  end=' ')
        for thing in self.corr:
            print(thing[0], end= ' ')
        print()


    def writeSpecTree(self):
        fh = open("../simData/spec.txt", "w")
        treeStr = self.sp_tree.as_string(schema='newick').split()[1]
        treeFinal = ''
        for c in treeStr:
            if c != "'":
                 treeFinal += c
        fh.write(treeFinal+"\n")
        fh.close()

    def runNullSims(self, pgls=False,RF=False,saveSumStats=False):
        
        rfs = []
        
        for k in range(self.numSims):
            gene_tree_new = treesim.contained_coalescent_tree(containing_tree=self.sp_tree,gene_to_containing_taxon_map=self.mapDict)

            rf0 = self.calcRF(gene_tree_new, self.sp_tree)
            rf1 = self.calcRF(self.gene_tree, gene_tree_new)

            if not RF:

                pdmGT = gene_tree_new.phylogenetic_distance_matrix()
                if  type(self.pdmST) == list:
                    self.pdmST = self.sp_tree.phylogenetic_distance_matrix()

                STdists = []
                GTdists = []
                trDist = []
                spec = []
         
                for node in self.gene_tree.leaf_node_iter():
                    if not pgls:
                        STdists = []
                        GTdists = []
                        trDist = []
                        spec = []

                    for node2 in self.gene_tree.leaf_node_iter():
                        if node2.taxon.label != node.taxon.label:
                            STdists.append( self.pdmST.distance(self.mapDict[node.taxon], self.mapDict[node2.taxon]))
                            GTdists.append( pdmGT.distance(self.mapDict[node.taxon], self.mapDict[node2.taxon]))
                            if self.trDist:
                                trDist.append(self.trDist[str(node.taxon)[1:-1]][str(node2.taxon)[1:-1]])                
                            else:
                                trDist.append((node.value-node2.value)**2)
                            spec.append([str(node.taxon),str(node2.taxon)])
                    res = linregress(GTdists, trDist)
                    if not pgls:
                        if not saveSumStats:
                            print(rf0, rf1, res.rvalue, res.pvalue, str(node.taxon)[1:-1])
                        else:
                            if str(node.taxon) in self.nullCorr:
                                self.nullCorr[str(node.taxon)].append(res.rvalue)
                            else:
                                self.nullCorr[str(node.taxon)] = [res.rvalue]
                if pgls:
                    for i in range(len(STdists)):
                        self.vals[(spec[i][0],spec[i][1])].append(GTdists[i])
            if RF:
               rfs.append(rf0) 

        if RF:
            return [np.mean(rfs), np.std(rfs)]


    def runNullSimsSet(self,RF=False):
        
        rfs = []
        
        null = []

        for k in range(self.numSims):
            gene_tree_new = treesim.contained_coalescent_tree(containing_tree=self.sp_tree,gene_to_containing_taxon_map=self.mapDict)

            rf0 = self.calcRF(gene_tree_new, self.sp_tree)
            rf1 = self.calcRF(self.gene_tree, gene_tree_new)

            if not RF:

                pdmGT = gene_tree_new.phylogenetic_distance_matrix()
                if  type(self.pdmST) == list:
                    self.pdmST = self.sp_tree.phylogenetic_distance_matrix()

                STdists = []
                GTdists = []
                trDist = []
                spec = []
         
                for node in self.gene_tree.leaf_node_iter():
                    for node2 in self.gene_tree.leaf_node_iter():
                        if node2.taxon.label != node.taxon.label:
                            STdists.append( self.pdmST.distance(self.mapDict[node.taxon], self.mapDict[node2.taxon]))
                            GTdists.append( pdmGT.distance(self.mapDict[node.taxon], self.mapDict[node2.taxon]))
                            if self.trDist:
                                trDist.append(self.trDist[str(node.taxon)[1:-1]][str(node2.taxon)[1:-1]])                
                            else:
                                trDist.append((node.value-node2.value)**2)
                            spec.append([str(node.taxon),str(node2.taxon)])
                res = linregress(GTdists, trDist)
                null.append(res.rvalue)

        nullMax = max(null)

        statPow = 0
        for val in self.corr:
            if val[0] >  nullMax:
                statPow += 1
        print(statPow,len(self.corr))


    def pglsPrint(self):
        k = -1
        old = ''
        for specPair in self.vals:
            new = specPair[0]
            if new != old:
                old = new
                k += 1
            print (str(specPair[0])[1:-1],str(specPair[1])[1:-1]+'-'+str(k),k,end=' ')
            for tr in self.vals[specPair]:
                print (tr, end= ' ')
            print()

    def pglsSims(self):
        
        self.setPopSize()
        self.makeTaxonMap()
        self.simGeneTree()
        self.simCausalGene(pgls=True)  
        self.runNullSims(pgls=True)       
        self.pglsPrint()        

    def simTraitAndNull(self):
        self.setPopSize()
        self.makeTaxonMap()
        self.simGeneTree()
        self.simCausalGene()
        self.writeSpecTree()
        self.makeOutputHeader()
        self.runNullSims()

    def simForRFDist(self,popSize):
        self.setPopSize(popSize=popSize)
        self.makeTaxonMap()
        self.simGeneTree()
        self.simCausalGene()
        self.writeSpecTree()
        myArr = self.runNullSims(RF=True)
        return myArr 

    def simRFDist(self,popSizeArr,empRF,xValsIn=[]):

        if not xValsIn:
            xVals = []
            for N in popSizeArr:
                xVals.append(self.simForRFDist(popSize=N)[0])
        else:
            xVals = xValsIn
   
        f = interp1d(xVals, -1.*np.log(popSizeArr),fill_value='extrapolate')
        retRF = []
        for r in empRF:
            minLogN = f(r)
            infN =  np.exp(-minLogN)
            retRF.append(infN)
        return retRF

    def setTrDist(self,trDist = {}):
        if not trDist:
            trDist = self.trDist 
        for node in self.gene_tree.leaf_node_iter():
             node.value = trDist[str(node.taxon)[1:-1]] 

    def getTrDist(self, f = '../obsData/RCCdata.txt'):
        fh = open(f, 'r')
        trDistLoc = defaultdict(dict)
        for line in fh:
            line = line.strip().split()
            if line[0] in trDistLoc and line[1] in trDistLoc[line[0]]: 
                trDistLoc[line[0]][line[1]].append(float(line[2]))
                trDistLoc[line[1]][line[0]].append(float(line[2]))
            else:
                trDistLoc[line[0]][line[1]] = [float(line[2])]
                trDistLoc[line[1]][line[0]] = [float(line[2])]

        # average the traits
        for s1 in trDistLoc:
            for s2 in trDistLoc[s1]:
                if np.mean(trDistLoc[s1][s2][0:3]) > np.mean(trDistLoc[s1][s2][3:]):
                    self.trDist[s1][s2] = np.mean(trDistLoc[s1][s2][0:3])  
                    self.trDist[s2][s1] = np.mean(trDistLoc[s1][s2][0:3])  
                else:           
                    self.trDist[s1][s2] = np.mean(trDistLoc[s1][s2][3:])
                    self.trDist[s2][s1] = np.mean(trDistLoc[s1][s2][3:])

    def correctPValDist(self,popSize,trDist={},save=True):
        self.setPopSize(popSize=popSize)   
        self.makeTaxonMap()
        self.simGeneTree()
        if len(trDist) == 0 and len(self.trDist) == 0:
            self.simCausalGene()
        else: 
            self.setTrDist(trDist)
        if save:
            self.runNullSims(saveSumStats=True)
        else:
            self.runNullSims(saveSumStats=False)

    def calcRealPVals(self):
        real_corr = []
        for tree in self.myTrees:
        
            pdmGT = tree.phylogenetic_distance_matrix()
            locCorr = []
            for node in tree.leaf_node_iter():
                GTdists = []
                trDist = []
                spec = []
                for node2 in tree.leaf_node_iter():
                    if str(node2.taxon) != str(node.taxon): # or node2.taxon == node.taxon:
                        GTdists.append( pdmGT.distance(node.taxon, node2.taxon))
                        trDist.append(self.trDist[str(node.taxon)[1:-1]][str(node2.taxon)[1:-1]])
                        spec.append([str(node.taxon)[1:-1],str(node2.taxon)[1:-1]])
                out = linregress(GTdists, trDist)
                locCorr.append( [out.rvalue, out.pvalue,str(node.taxon)[1:-1]])
            real_corr.append(locCorr) 
        if real_corr:
            self.real_corr = real_corr

    def readGeneTree(self,myDir="/cluster/tufts/uricchiolab/nlouw01/post_gensas/OrthoFinder/Results_Mar14/Gene_Trees",fn=1):
        zers = 7-len(str(fn))
        zrStr = ''
        for i in range(zers):
            zrStr += '0'
        tree = dendropy.Tree.get(path=os.path.join(myDir,'OG'+zrStr+str(fn)+'_tree.txt'), schema="newick",taxon_namespace=self.sp_tree.taxon_namespace)
        return tree

    def getSingleCopy(self,myMin=1,myMax = 9349, specList = set(['Penicillium-polonicum','Penicillium-verrucosum','Penicillium-solitum','Penicillium-mb',
                                                     'Penicillium-cvjetkovicii','Penicillium-commune','Penicillium-chrysogenum','Penicillium-brevicompactum'])):
        copyNum = []
        for i in range(myMin,myMax+1):
            zers = 7-len(str(i))
            zerStr = ''
            for j in range(zers):
                zerStr += '0'
            fasta = "/cluster/tufts/uricchiolab/nlouw01/post_gensas/OrthoFinder/Results_Mar14/Single_Copy_Orthologue_Sequences/OG"+zerStr+str(i)+".fa"
            if not os.path.isfile(fasta):
                continue
            myTree = self.readGeneTree(fn=i)
            k = 0
            names = set()
            for l in myTree.leaf_node_iter():
                name = str(l.taxon).split()[0][1:-3]
                l.label = name
                l.taxon = taxonmodel.Taxon(name)
                names.add(name)
                k += 1
                if k > len(specList):
                    break
            #if names == specList and k <= len(specList):
            myTree.migrate_taxon_namespace(self.sp_tree.taxon_namespace)
            self.myTrees.append(myTree)
            copyNum.append(i)
        self.copyNum = copyNum  

    def calcRFDist(self):
    
        rfs = []
        for tr in self.myTrees:
            tr.encode_bipartitions()
            
            rfs.append(treecompare.symmetric_difference(tr, self.sp_tree))
        return [np.mean(rfs), np.std(rfs)] 

    def getIndVals(self,geneNum=5719):    
        i = 0
        for tree in self.myTrees:
            if self.copyNum[i] != geneNum:
                i += 1
                continue 
  
            pdmGT = tree.phylogenetic_distance_matrix()
            pdmST = self.sp_tree.phylogenetic_distance_matrix()
            locCorr = []
            for node in tree.leaf_node_iter():
                GTdists = []
                STdists = []
                trDist = []
                spec = []
                for node2 in tree.leaf_node_iter():
                    if str(node2.taxon) != str(node.taxon): # or node2.taxon == node.taxon:
                        GTdists.append( pdmGT.distance(node.taxon, node2.taxon))
                        STdists.append( pdmST.distance(node.taxon, node2.taxon))
                        trDist.append(self.trDist[str(node.taxon)[1:-1]][str(node2.taxon)[1:-1]])
                        spec.append([str(node.taxon)[1:-1],str(node2.taxon)[1:-1]])
                
                out = linregress(GTdists, trDist)
                print("#", out.rvalue, out.pvalue,str(node.taxon)[1:-1])
                for j in range(len(GTdists)):
                    print(GTdists[j], trDist[j], str(node.taxon)[1:-1], STdists[j])
            break          
      
    def getCorrPVals(self,nullDistFile='../obsData/corrPvalDist.txt',trueCorrCoeffFile='../obsData/realPvalDist.txt'):
   
        nullDist = {}   
 
        # get null dist for each species
        fh = open(nullDistFile, 'r')
        for line in fh:
            data = line.strip().split()
            if data[0] not in nullDist:
                nullDist[data[0]] = []
            nullDist[data[0]].append(float(data[1]))

        fh.close() 

        dist = defaultdict(dict)
        # for each gene tree, get real corr coeff
        fh = open(trueCorrCoeffFile, 'r')
        for line in fh:
            data = line.strip().split()
            dist[data[2]][data[3]] = float(data[0])       
    
        # for each corr coeff, compute quantile
        percentiles = {}
        gtOrder = {}
        for spec in nullDist:
            nullDist[spec] = sorted(nullDist[spec])
            allCC = []
            for gt in dist[spec]:
                allCC.append([dist[spec][gt],gt])
            allCC = sorted(allCC, key=lambda x: x[0])
            i = 0
            percentiles[spec] = []
            
            gtOrder[spec] = []
            while allCC:
                if len(nullDist[spec]) > i and nullDist[spec][i] < allCC[0][0]:
                    i += 1
                    continue
                elif len(nullDist[spec]) <= i:
                    for item in allCC:
                        percentiles[spec].append(0.)
                        gtOrder[spec].append(item[1])
                    break    
                val = allCC.pop(0)
                percentiles[spec].append(1.-i/(len(nullDist[spec])))                 
                gtOrder[spec].append(val[1])
 
        return (percentiles,gtOrder)              

    def simCausalGeneSet(self):
        corr = []
        for gt in self.gene_tree_list:
            pdmGT = gt.phylogenetic_distance_matrix()
            self.pdmST = self.sp_tree.phylogenetic_distance_matrix()

            STdists = []
            GTdists = []
            trDist = []
            for node in gt.leaf_node_iter():
                spec = []
                STdists = []
                for node2 in gt.leaf_node_iter():
                    if node2.taxon.label != node.taxon.label: # or node2.taxon == node.taxon:
                        STdists.append( self.pdmST.distance(self.mapDict[node.taxon], self.mapDict[node2.taxon]))
                        GTdists.append( pdmGT.distance(self.mapDict[node.taxon], self.mapDict[node2.taxon]))
                        trDist.append((node.value-node2.value)**2)
                        spec.append([str(node.taxon),str(node2.taxon)])
            out = linregress(GTdists, trDist)
            corr.append( [out.rvalue, out.pvalue])
        self.corr = corr

    def simPowerSpecTree(self,nLin=8,nGene=1):
        
        # first, simulate species tree with Birth Death process
        self.initRandomTree(n=nLin)
        #print(self.sp_tree.as_string(schema="newick"))

        # now simualte underlying causal genes
        self.setPopSize()
        self.makeTaxonMap()
        for i in range(nGene):
            self.simGeneTree()
            self.evolveTrait()
            self.gene_tree = self.normalize(self.gene_tree)
            self.gene_tree_list.append(self.gene_tree) 
        trDist = {} 
        for tr in self.gene_tree_list:
            for node in tr.leaf_node_iter():
                if node.taxon.label in trDist:
                     trDist[node.taxon.label] += node.value
                else:
                     trDist[node.taxon.label]= node.value
        for node in self.sp_tree.leaf_node_iter():
            node.value = trDist[node.taxon.label]

        for tr in self.gene_tree_list:
            for node in tr.leaf_node_iter():
                node.value = trDist[node.taxon.label]

        self.simCausalGeneSet()
        self.runNullSimsSet()

        return

    def calcSumStats(self,param,simDistFile,realDistFile):
         
        simDist = {}
        realDist = {} 

        #read simDist
        fh = open(simDistFile,'r')
        for line in fh:
            data = line.strip().split() 
            if data[0] in simDist:      
                simDist[data[0]].append(float(data[1]))
            else:
                simDist[data[0]] = [float(data[1])]
        fh.close()

        # read realDist
        fh = open(realDistFile,'r')
        for line in fh:
            data = line.strip().split()
            if data[2] in realDist:
                realDist[data[2]].append(float(data[0]))
            else:
                realDist[data[2]] = [float(data[0])]
        fh.close()

        # for spec in simDist, calc diff in mean and std and sum them
        tot = 0
        for spec in realDist:
            mr = np.mean(realDist[spec])
            sr = np.std(realDist[spec])

            ms = np.mean(simDist[spec])        
            ss = np.std(simDist[spec]) 

            tot += (ms-mr)**2 + (ss-sr)**2
                
        # print param, distance
        print(param, tot)

    def interpPopSize(self,popSizeStatFile):
        
        xvals = []
        yvals = []
        fh = open(popSizeStatFile, 'r')
        for line in fh:
            data = line.strip().split()
            xvals.append(float(data[0]))
            yvals.append(float(data[1]))
        fh.close()

        f = interp1d(xvals, yvals,fill_value='extrapolate')
        res  = minimize(f,x0=0.083).x[0]

        print(res)

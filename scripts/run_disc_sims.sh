# ILS/correlation analysis
#for i in {1..1000}; do python sim.py 100 100000 > ../simData/highDiscP.${i}.txt; done
#for i in {1..1000}; do python sim.py 100 0.001 > ../simData/lowDiscP.${i}.txt; done
#for i in {1..1000}; do python sim.py 100 0.1 > ../simData/midDiscP.${i}.txt; done

# ILS/correlation analysis group of genes
#./subBatchAll.sh

# PGLS analysis
#for i in {1..100}; do python simPGLS.py 100 100000 > ../PGLSsims/pgls.highDiscP.${i}.txt; done
#for i in {1..100}; do python simPGLS.py 100 0.001 > ../PGLSsims/pgls.lowDiscP.${i}.txt; done
#for i in {1..100}; do python simPGLS.py 100 0.1 > ../PGLSsims/pgls.midDiscP.${i}.txt; done

# calc corr pVals with simulated trees
#for pS in {0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12}; do python pValsNull.py $pS 10000  > ../obsData/pValsNull/pValsNull.$pS.txt; done

# get Sum Stat diffs for inferring N
#for i in {0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12}; do python getSumStatDiff.py $i >> ../obsData/pValsNull/sumDataDiffs.txt; done

# minimize to get popSize for simulated null
python minimize.py > ../obsData/popSize.txt

#get Null pVals for correct pop Size
python pVals.py 1000000  > ../obsData/corrPvalDist.txt

# calc corr pVals, not saving sum stats
python pValsCorrBtwn.py 100000  > ../obsData/corrPvalDistCorrBtwn.txt

#calc dist of real pVals and corr Coeffs
python pValsReal.py > ../obsData/realPvalDist.txt

# data for 5719
python getDataForVis.py > ../obsData/corr.5719.txt

# get list of genes sig at p < alpha level with X+ assoc, results for table S3
#cat ../obsData/QQplotData.txt | awk '{if($4 < 0.05) print $3}' | sort | uniq -c | awk '{if($1 > 2) print $2}' | wc -l
#cat ../obsData/QQplotData.txt | awk '{if($4 < 0.05) print $3}' | sort | uniq -c | awk '{if($1 > 3) print $2}' | wc -l
#cat ../obsData/QQplotData.txt | awk '{if($4 < 0.01) print $3}' | sort | uniq -c | awk '{if($1 > 2) print $2}' | wc -l
#cat ../obsData/QQplotData.txt | awk '{if($4 < 0.01) print $3}' | sort | uniq -c | awk '{if($1 > 3) print $2}' | wc -l

#results of sims for Table S3
#python getQuantiles.py ../obsData/corrPvalDistCorrBtwn.txt 0.01 3 
#python getQuantiles.py ../obsData/corrPvalDistCorrBtwn.txt 0.01 4 
#python getQuantiles.py ../obsData/corrPvalDistCorrBtwn.txt 0.05 3 
#python getQuantiles.py ../obsData/corrPvalDistCorrBtwn.txt 0.05 4 

#making enrichment gene sets for kegg analysis
#cat ../../../QQplotData.txt | awk '{if($4 < 0.005) print $3}' | sort | uniq -c | awk '{if($1 > 0) print $2}' > ../obsData/compHits.plt0.005.1sp.txt 
#cat ../../../QQplotData.txt | awk '{if($4 < 0.0005) print $3}' | sort | uniq -c | awk '{if($1 > 0) print $2}' > ../obsData/compHits.plt0.0005.1sp.txt 
#cat ../../../QQplotData.txt | awk '{if($4 < 0.05) print $3}' | sort | uniq -c | awk '{if($1 > 2) print $2}' > ../obsData/compHits.plt0.05.gt2sp.txt 
#cat ../../../QQplotData.txt | awk '{if($4 < 0.05) print $3}' | sort | uniq -c | awk '{if($1 > 3) print $2}' > ../obsData/compHits.plt0.05.gt3sp.txt 
#cat ../../../QQplotData.txt | awk '{if($4 < 0.01) print $3}' | sort | uniq -c | awk '{if($1 > 2) print $2}' > ../obsData/compHits.plt0.01.gt2sp.txt 
#cat ../../../QQplotData.txt | awk '{if($4 < 0.01) print $3}' | sort | uniq -c | awk '{if($1 > 3) print $2}' > ../obsData/compHits.plt0.01.gt3sp.txt 

#python getIDs.py ../obsData/GeneMap.txt ../obsData/compHits.plt0.005.1sp.txt > ../obsData/enrichment/compHits.plt0.005.1sp.txt
#python getIDs.py ../obsData/GeneMap.txt ../obsData/compHits.plt0.0005.1sp.txt > ../obsData/enrichment/compHits.plt0.0005.1sp.txt
#python getIDs.py ../obsData/GeneMap.txt ../obsData/compHits.plt0.05.gt2sp.txt > ../obsData/enrichment/compHits.plt0.05.gt2sp.txt
#python getIDs.py ../obsData/GeneMap.txt ../obsData/compHits.plt0.05.gt3sp.txt > ../obsData/enrichment/compHits.plt0.05.gt3sp.txt
#python getIDs.py ../obsData/GeneMap.txt ../obsData/compHits.plt0.01.gt2sp.txt > ../obsData/enrichment/compHits.plt0.01.gt2sp.txt
#python getIDs.py ../obsData/GeneMap.txt ../obsData/compHits.plt0.01.gt3sp.txt > ../obsData/enrichment/compHits.plt0.01.gt3sp.txt


# not using any of this below anymore
# Binary trait regression
#for i in {1..200}; do python binTraitReg.py 100 100000 > ../phyloGWASsims/bin.highDiscP.${i}.txt; done
#for i in {1..200}; do python binTraitReg.py 100 0.001 > ../phyloGWASsims/bin.lowDiscP.${i}.txt; done
#for i in {1..200}; do python binTraitReg.py 100 0.1 > ../phyloGWASsims/bin.midDiscP.${i}.txt; done

# phyloGWAS
#for i in {101..250}; do python phyloGWASNull.py 1000 100000 > ../phyloGWASsims/pGWAS.highDiscP.${i}.txt; done
#for i in {101..250}; do python phyloGWASNull.py 1000 0.001 > ../phyloGWASsims/pGWAS.lowDiscP.${i}.txt; done
#for i in {101..250}; do python phyloGWASNull.py 1000 0.1 > ../phyloGWASsims/pGWAS.midDiscP.${i}.txt; done

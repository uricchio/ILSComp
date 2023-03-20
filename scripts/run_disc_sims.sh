# ILS/correlation analysis
#for i in {1..1000}; do python sim.py 100 100000 > ../simData/highDiscP.${i}.txt; done
#for i in {1..1000}; do python sim.py 100 0.001 > ../simData/lowDiscP.${i}.txt; done
#for i in {1..1000}; do python sim.py 100 0.1 > ../simData/midDiscP.${i}.txt; done

# PGLS analysis
#for i in {1..100}; do python simPGLS.py 100 100000 > ../PGLSsims/pgls.highDiscP.${i}.txt; done
#for i in {1..100}; do python simPGLS.py 100 0.001 > ../PGLSsims/pgls.lowDiscP.${i}.txt; done
#for i in {1..100}; do python simPGLS.py 100 0.1 > ../PGLSsims/pgls.midDiscP.${i}.txt; done

# inferring pop size for null dist of corr
#python infer.py 10000 > ../simData/popSizeInfer.txt

# calc dist of RF on real trees
#python calcTreeDisc.py > ../obsData/rfDistGT.txt

# calc corr pVals with simulated trees
#python pVals.py 10000 1000000  > ../obsData/corrPvalDist.txt

# calc corr pVals, not saving sum stats
#python pValsCorrBtwn.py 10000 100000  > ../obsData/corrPvalDistCorrBtwn.txt

#calc dist of real pVals and corr Coeffs
#python pValsReal.py > ../obsData/realPvalDist.txt

# data for 5719
#python getDataForVis.py > ../obsData/corr.5719.txt

# not using any of this below anymore
# Binary trait regression
#for i in {1..200}; do python binTraitReg.py 100 100000 > ../phyloGWASsims/bin.highDiscP.${i}.txt; done
#for i in {1..200}; do python binTraitReg.py 100 0.001 > ../phyloGWASsims/bin.lowDiscP.${i}.txt; done
#for i in {1..200}; do python binTraitReg.py 100 0.1 > ../phyloGWASsims/bin.midDiscP.${i}.txt; done

# phyloGWAS
#for i in {101..250}; do python phyloGWASNull.py 1000 100000 > ../phyloGWASsims/pGWAS.highDiscP.${i}.txt; done
#for i in {101..250}; do python phyloGWASNull.py 1000 0.001 > ../phyloGWASsims/pGWAS.lowDiscP.${i}.txt; done
#for i in {101..250}; do python phyloGWASNull.py 1000 0.1 > ../phyloGWASsims/pGWAS.midDiscP.${i}.txt; done

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

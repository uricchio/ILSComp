#cat ../obsData/delta_anti_merged2.csv | awk '{FS=","; print $2,$3,$4}' | awk '{if ($1 != "Alone" && $2 != "Alone") print}' | awk '{if($3 < 0) print $1,$2,-1*$3; else print}' > ../obsData/delta_anti_merged2.txt

# get null dist of corr coeffs
#python pValsRuby.py 10000 1000000  > ../obsData/corrPvalDistRuby.txt

# get readl corr coeffs
#python pValsRealRuby.py > ../obsData/realPvalDistRuby.txt

# get est p-values under null
#python getQuantilesReal.py > ../obsData/correctedPvalsRuby.txt

sed 's/^/chr/' dominantGenesRegions.bed | awk '$4=NR {print}' > dominantGenesRegions_chr.bed
sed 's/^/chr/' non.dominantGenesRegions.bed | awk '$4=NR {print}' > non.dominantGenesRegions_chr.bed
mafsInRegion dominantGenesRegions_chr.bed -outDir dominant *.maf
mafsInRegion non.dominantGenesRegions_chr.bed -outDir nondominant *.maf

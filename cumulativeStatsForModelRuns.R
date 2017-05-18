cumBiomass = c(
  CN_anoxic = sum(sapply(1:9, function(i) resultsCN_Ox.Hx [[i]][[1]])),
  CN_oxic = sum(sapply(1:9, function(i) resultsCN_Ox.O2.Hx [[i]][[1]])),
  CNO = sum(sapply(1:9, function(i) resultsCON_Ox.Hx [[i]][[1]])),
  CNOS = sum(sapply(1:9, function(i) resultsCONSH_Ox.Hx[[i]][[1]]))
)


cumEnergies = c(
  CN_anoxic = sum(sapply(1:9, function(i) sum(resultsCN_Ox.Hx[[i]]$processEnergyVals$energy[resultsCN_Ox.Hx[[i]]$processEnergyVals$procType == "catabolic"]))),
  CN_oxic = sum(sapply(1:9, function(i) sum(resultsCN_Ox.O2.Hx[[i]]$processEnergyVals$energy[resultsCN_Ox.O2.Hx[[i]]$processEnergyVals$procType == "catabolic"]))),
  CNO = sum(sapply(1:9, function(i) sum(resultsCON_Ox.Hx[[i]]$processEnergyVals$energy[resultsCON_Ox.Hx[[i]]$processEnergyVals$procType == "catabolic"]))),
  CNOS = sum(sapply(1:9, function(i) sum(resultsCONSH_Ox.Hx[[i]]$processEnergyVals$energy[resultsCONSH_Ox.Hx[[i]]$processEnergyVals$procType == "catabolic"])))
)


par(mfrow = c(2, 1))
barplot(cumEnergies)
barplot(cumBiomass)

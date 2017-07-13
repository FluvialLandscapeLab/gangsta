cumBiomass = c(
  CN_anoxic = mean(sapply(1:9, function(i) resultsCN_Ox.Hx [[i]][[1]])),
  CN_oxic = mean(sapply(1:9, function(i) resultsCN_Ox.O2.Hx [[i]][[1]])),
  CNO = mean(sapply(1:9, function(i) resultsCON_Ox.Hx [[i]][[1]])),
  CNOS = mean(sapply(1:9, function(i) resultsCONSH_Ox.Hx[[i]][[1]]))
)


cumEnergies = c(
  CN_anoxic = mean(sapply(1:9, function(i) sum(resultsCN_Ox.Hx[[i]]$processEnergyVals$energy[resultsCN_Ox.Hx[[i]]$processEnergyVals$procType == "catabolic"]))),
  CN_oxic = mean(sapply(1:9, function(i) sum(resultsCN_Ox.O2.Hx[[i]]$processEnergyVals$energy[resultsCN_Ox.O2.Hx[[i]]$processEnergyVals$procType == "catabolic"]))),
  CNO = mean(sapply(1:9, function(i) sum(resultsCON_Ox.Hx[[i]]$processEnergyVals$energy[resultsCON_Ox.Hx[[i]]$processEnergyVals$procType == "catabolic"]))),
  CNOS = mean(sapply(1:9, function(i) sum(resultsCONSH_Ox.Hx[[i]]$processEnergyVals$energy[resultsCONSH_Ox.Hx[[i]]$processEnergyVals$procType == "catabolic"])))
)


par(
  mfrow = c(2, 1),
  mar = c(5.5, 5, 4, 2) + 0.1
  )
barplot(
  cumBiomass,
  ylab = "", xlab = "", names.arg = c("CN Anoxic", "CN Oxic", "CNO", "CNOS"), main= "Biomass",
  cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2,
  las = 1
  )
title(ylab = expression(paste("\u03bc", "mol")), line = 3.75, cex.lab = 1.2)
barplot(
  cumEnergies*1000,
  ylab = "", xlab = "", names.arg = c("CN Anoxic", "CN Oxic", "CNO", "CNOS"), main = "Dissimilatory Energy",
  cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2,
  las = 1
  )
title(ylab = "Joules", line = 3.75, cex.lab = 1.2)



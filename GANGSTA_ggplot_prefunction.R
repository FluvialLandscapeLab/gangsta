# source('M:/gangsta/R/GANGSTA_ggplot_input_centering_addToPackage.R')


CNOSH_Any(c("C", "N"),
          c("Ox", "Hx"))
CNOSH_Any(c("C", "N"),
          c("Ox", "O2", "Hx"))
CNOSH_Any(c("C", "O", "N"),
          c("Ox", "Hx"))
CNOSH_Any(c("C", "O", "N", "S", "H"),
          c("Ox", "Hx"))

makePlots("resultsCON_Ox.Hx", "gangstasCON_Ox.Hx", elementalCyclesToPlot = "O", aggregateBio = F, yAxisMaxMols = 0, axisFontSize = 1.5)
makePlots("resultsCONSH_Ox.Hx", "gangstasCONSH_Ox.Hx", elementalCyclesToPlot = "O", aggregateBio = F, yAxisMaxMols = 0, axisFontSize = 1.5)

# makeOutput()
makePlots(pdf = F, aggregateBio = F, yAxisMaxMols = 0, axisFontSize = 1.5)
# makePlots(pdf = T, aggregateBio = F, elementalCyclesToPlot ="C")
# makePlots(pdf = T, aggregateBio = F, elementalCyclesToPlot ="N")
# makePlots(pdf = F, aggregateBio = F, elementalCyclesToPlot = c("C", "N"))

combineRiverPlotsInPDF(
  fileIdx = "CNOS_TEST_20170301",
  cyclesSeparate = T,
  elementalCyclesToPlot = c("C", "N", "O", "S"),
  axisFontSize = 0.5,
  yAxisMaxMols = 4.05
  )



# ############### DISSIM ENERGY PLOTS


combineDissimEnergyPlotsInPDF(withLegend = F, fileIdx = "dissimPlots_noLegend_20170424", axisFontSize = 1.5)
combineDissimEnergyPlotsInPDF(withLegend = T, fileIdx = "dissimPlots_withLegend_20170410")


########## Energy balance plots
combineEnergyBalPlotsInPDF(withLegend = F, fileIdx = "energyBalPlots_noLegend_20170424", axisFontSize = 1.5)


#### Biomass plots
combineBiomassPlotsInPDF(withLegend = F, fileIdx = "biomassPlots_noLegend_20170424", axisFontSize = 1.5)


### Substrates and products plots
combineSubstrProdPlots(fileIdx = "reactantsAndProductsPlots_noLegend20170424", axisFontSize = 1.5)

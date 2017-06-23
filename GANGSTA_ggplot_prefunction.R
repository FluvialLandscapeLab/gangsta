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
  fileIdx = "CNOS_TEST_20170623",
  cyclesSeparate = T,
  elementalCyclesToPlot = c("C", "N", "O", "S"),
  axisFontSize = 0.75,
  yAxisMaxMols = 4.05
  )

########## Energy balance plots
combineEnergyBalPlotsInPDF(withLegend = F, fileIdx = "energyBalPlots_noLegend_20170623", axisFontSize = 1.2)


#### Biomass plots
combineBiomassPlotsInPDF(withLegend = F, fileIdx = "biomassPlots_noLegend_20170623", axisFontSize = 1.2)

# ############### DISSIM ENERGY PLOTS
combineDissimEnergyPlotsInPDF(withLegend = F, fileIdx = "dissimPlots_noLegend_20170623", axisFontSize = 1.2)

### combine leakIn  plots
combineLeakInPlots(fileIdx = "leakInPlots_noLegend20170620", axisFontSize = 1.2)
combineCompleteLeakInListPlots(fileIdx = "completeLeakInPlots_noLegend20170623", axisFontSize = 1.2)

### combine product  plots
combineProductsPlots(fileIdx = "ProductsPlots_noLegend20170623", axisFontSize = 1.2)

# I used the following to get legends of the same size for the dissim and leak in/prod plots:
makeDissimEnergyPlot(resultsList = resultsCONSH_Ox.Hx, withLegend = T, printPDF = T, textCol = "black", backgroundCol = "white", axisFontSize = 1.2)
productsPlot(resultsList = resultsCONSH_Ox.Hx, withLegend = T, printPDF = T, textCol = "black", backgroundCol = "white", axisFontSize = 1.2)

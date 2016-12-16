# source('M:/gangsta/R/GANGSTA_ggplot_input_centering_addToPackage.R')


CNOSH_Any(c("C", "N"),
          c("Ox", "Hx"))
CNOSH_Any(c("C", "N"),
          c("Ox", "O2", "Hx"))
CNOSH_Any(c("C", "O", "N"),
          c("Ox", "Hx"))
CNOSH_Any(c("C", "O", "N", "S", "H"),
          c("Ox", "Hx"))


# makeOutput()
makePlots(pdf = T, aggregateBio = F, yAxisMaxMols = 0, axisFontSize = 1.5)
# makePlots(pdf = T, aggregateBio = F, elementalCyclesToPlot ="C")
# makePlots(pdf = T, aggregateBio = F, elementalCyclesToPlot ="N")
# makePlots(pdf = F, aggregateBio = F, elementalCyclesToPlot = c("C", "N"))

combineRiverPlotsInPDF(
  fileIdx = "CNOS_TEST_20161212",
  cyclesSeparate = T,
  elementalCyclesToPlot = c("C", "N", "O", "S"),
  axisFontSize = 0.5,
  yAxisMaxMols = 4.05
  )



# ############### DISSIM ENERGY PLOTS


combineDissimEnergyPlotsInPDF(withLegend = F, fileIdx = "dissimPlots_noLegend_12122016", axisFontSize = 2)
combineDissimEnergyPlotsInPDF(withLegend = T, fileIdx = "dissimPlots_withLegend_11292016")

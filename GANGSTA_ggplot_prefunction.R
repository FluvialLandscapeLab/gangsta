source('M:/gangsta/GANGSTA_ggplot_input_centering.R')
source('M:/gangsta/leakInList.R')




CNOSH_Any(c("C", "N"),
          c("Ox", "Hx"))
CNOSH_Any(c("C", "N"),
          c("Ox", "O2", "Hx"))
CNOSH_Any(c("C", "O", "N"),
          c("Ox", "Hx"))
CNOSH_Any(c("C", "O", "N", "S", "H"),
          c("Ox", "Hx"))


# makeOutput()
# makePlots(pdf = T, aggregateBio = F)
# makePlots(pdf = T, aggregateBio = F, elementalCyclesToPlot ="C")
# makePlots(pdf = T, aggregateBio = F, elementalCyclesToPlot ="N")
# makePlots(pdf = F, aggregateBio = F, elementalCyclesToPlot = c("C", "N"))

combineRiverPlotsInPDF(
  fileIdx = "CNOS_CombinedPlots_sep",
  cyclesSeparate = T,
  elementalCyclesToPlot = c("C", "N", "O", "S"),
  axisFontSize = 0.5
  )



# ############### DISSIM ENERGY PLOTS


combineDissimEnergyPlotsInPDF(withLegend = F, fileIdx = "dissimPlots_noLegend")
combineDissimEnergyPlotsInPDF(withLegend = T, fileIdx = "dissimPlots_withLegend")



#
# ############### DISSIM ENERGY PLOTS
# combineDissimEnergyPlotsInPDF = function(withLegend = F, fileIdx){
#
#   filePrefix = "C:\\Users\\AnnMarie\\Documents\\Research\\Projects\\BGC\\SFS 2016\\DissimEnergyPlots\\"
#   fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)
#   dissimPlot.C =
#     makeDissimEnergyPlot(results = resultsC_Ox.Hx, withLegend = withLegend)
#   dissimPlot.C.aerobic =
#     makeDissimEnergyPlot(results = resultsC_O2.Ox.Hx, withLegend = withLegend)

#
#   dissimPlotList =
#     list(
#       dissimPlot.C,
#       dissimPlot.C.aerobic,
#       dissimPlot.CO,
#       dissimPlot.CON,
#       dissimPlot.CHONS
#     )
#   pdf(
#     fileName,
#     onefile = T,
#     paper = "USr",
#     # paper = "letter"
#     width = 10.5,
#     height = 7.5
#   )
#   print(grid.arrange(grobs = dissimPlotList, ncol = 1))
#   dev.off()
# }
# combineDissimEnergyPlotsInPDF(withLegend = F, fileIdx = "dissimPlots_noLegend")
# combineDissimEnergyPlotsInPDF(withLegend = T, fileIdx = "dissimPlots_withLegend")
#
#










# plotList_superPlot.CO =
#   list(
#     superPlot.CO_C,
#     superPlot.CO_O,
#     superPlot.CO
#   )
# pdf(
#   "C:\\Users\\AnnMarie\\Documents\\Research\\Projects\\BGC\\SFS 2016\\SuperPlots\\superPlot.CO.pdf",
#   onefile = T,
#   # paper = "USr",
#   paper = "letter",
#   width = 8,
#    height = 10.5
# )
# print(grid.arrange(grobs = plotList_superPlot.CO, ncol = 1))
# dev.off()
#
#
# ##
#
# inputDF_test = inputDF.CHONS[(inputDF.CHONS$element %in% c("C", "O", "N", "S")),]
# gangstaSuperPlot(inputDF_test, makePDF = F)
#
#
#

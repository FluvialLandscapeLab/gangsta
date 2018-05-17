

CNOSH_Any(c("C", "N"),
          c("Ox", "Hx"))
CNOSH_Any(c("C", "N"),
          c("Ox", "O2", "Hx"))
CNOSH_Any(c("C", "O", "N"),
          c("Ox", "Hx"))
CNOSH_Any(c("C", "O", "N", "S", "H"),
          c("Ox", "Hx"))


### River plots
makePlots(
  pdf = T,
  "resultsCONSH_Ox.Hx",
  "gangstasCONSH_Ox.Hx",
  elementalCyclesToPlot = c("C", "H", "O", "N", "S"),
  aggregateBio = F,
  yAxisMaxMols = 0,
  axisFontSize = 1.5,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\RiverPlots\\")

combineRiverPlotsInPDF(
  fileIdx = "CNOS_TEST_20180118",
  cyclesSeparate = T,
  elementalCyclesToPlot = c("C", "N", "O", "S"),
  axisFontSize = 0.75,
  yAxisMaxMols = 4.05,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\RiverPlots\\"
  )



########## Energy balance plots
combineEnergyBalPlotsInPDF(
  withLegend = F,
  fileIdx = "energyBalPlots_noLegend_20180108",
  axisFontSize = 1.2,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\EnergyBalancePlots\\"
)


#### Biomass plots
combineBiomassPlotsInPDF(
  withLegend = F,
  fileIdx = "biomassPlots_noLegend_20180104",
  axisFontSize = 1.2,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\BiomassPlots\\"
)

# ############### DISSIM ENERGY PLOTS
combineDissimEnergyPlotsInPDF(
  withLegend = F,
  fileIdx = "dissimPlots_noLegend_20170108",
  axisFontSize = 1.2,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\DissimEnergyPlots\\"
)


### combine leakIn  plots
combineLeakInPlots(
  fileIdx = "leakInPlots_noLegend20180104",
  axisFontSize = 1.2,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\")
combineCompleteLeakInListPlots(
  fileIdx = "completeLeakInPlots_noLegend20180108",
  axisFontSize = 1.2,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\")

### combine product  plots
combineProductsPlots(
  fileIdx = "ProductsPlots_noLegend20180104",
  axisFontSize = 1.2,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\"
)

### combine substrate and product plots (Inputs, Reactants, and Products plots)
combineSubstrProdPlots(
  fileIdx = "SubstratesAndProductsPlots_TEST_20180104",
  axisFontSize = 1.2,
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\"
)

# I used the following to get legends of the same size for the dissim and leak in/prod plots:
makeDissimEnergyPlot(
  resultsList = resultsCONSH_Ox.Hx,
  withLegend = T,
  printPDF = T,
  textCol = "black",
  backgroundCol = "white",
  axisFontSize = 1.2,
  fileName = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\DissimEnergyPlots\\dissimEnergyPlot_20180104.pdf"
  )
productsPlot(
  resultsList = resultsCONSH_Ox.Hx,
  withLegend = T,
  printPDF = T,
  textCol = "black",
  backgroundCol = "white",
  axisFontSize = 1.2,
  fileName = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\prodPlot_20180104.pdf"
)

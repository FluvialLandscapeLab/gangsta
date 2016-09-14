makeDissimEnergyPlot = function(
  resultsList,
  printPDF = F,
  fileName = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\DissimEnergyPlots\\dissimEnergyPlot.pdf",
  withLegend = F,
  textCol,
  backgroundCol,
  axisFontSize
) {
  if(withLegend ==F) {
    legendLocation = "none"
  }else{
    legendLocation = "top"
  }

  dissimEnergyList = lapply(resultsList, function(x) x$processEnergyVals[x$processEnergyVals$procType == "catabolic", ])
  dissimEnergyDF = dissimEnergyList[[1]]
  dissimEnergyDF = data.frame(energy = dissimEnergyDF$energy, timestep = 1, process = row.names(dissimEnergyDF))
  for(i in 2:length(dissimEnergyList)){
    newRow = data.frame(energy = dissimEnergyList[[i]]$energy, timestep = i, process = row.names(dissimEnergyList[[i]]))
    dissimEnergyDF = rbind(dissimEnergyDF, newRow)
  }
  dissimEnergyDF$energyInJoules = dissimEnergyDF$energy * 1000

  processNames =
    c("AutNitrif", "AutSulfideOxidation", "HetAerobic", "HetDenit", "HetMethanogenesis", "HetSulfateRed", "MetMethaneOxid")
  processColors =
    c("chartreuse3", "darkorange", "blue3", "brown3", "purple", "yellow", "cyan3")
  processLabels =
    c("Nitrification", "Sulfide oxidation", "Aerobic heterotrophy", "Denitrification", "Methanogenesis", "Sulfate reduction", "Methane oxidation")
  names(processColors) = processNames
  names(processLabels) = processNames

  processColors = processColors[names(processColors) %in% levels(dissimEnergyDF$process)]
  processLabels = processLabels[names(processLabels) %in% levels(dissimEnergyDF$process)]

  dissimEnergyPlot = ggplot(dissimEnergyDF, aes(x = timestep, y = energyInJoules, fill = process)) +
    geom_bar(stat = "identity") +
    expand_limits(y=c(0,0.1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.background = element_rect(fill = backgroundCol, colour = backgroundCol),
          panel.background = element_blank()) +
    scale_x_continuous(
      breaks = seq(1,9,1),
      labels = as.character(seq(1, 9, 1)),
      name = element_blank()
    ) +
    scale_y_continuous(
      breaks = seq(0, 0.09, 0.03),
      labels = as.character(seq(0, 0.09, 0.03)),
      name = ""
    ) +
    theme(axis.text = element_text(colour = textCol, size = rel(axisFontSize))) +
    theme(
      legend.background = element_rect(fill = backgroundCol),
      legend.text = element_text(size = rel(1.2), colour = textCol  ),
      legend.position = legendLocation) +
    scale_fill_manual(
      values = processColors,
      # c("chartreuse3", "darkorange", "blue3", "brown3", "purple", "yellow", "cyan3"),
      labels = processLabels
      # c("Nitrification", "Sulfide oxidation", "Aerobic heterotrophy", "Denitrification", "Methanogenesis", "Sulfate reduction", "Methane oxidation")
    ) +
    theme(plot.margin = unit(c(-0.60, 1, 0.02, 0.5), "lines"))
  if(printPDF == TRUE) {
    pdf(
      fileName,
      onefile = T,
      paper = "USr",
      # paper = "letter",
      width = 11,
      height = 4.5
    )
    print(dissimEnergyPlot)
    dev.off()
  } else {
    return(dissimEnergyPlot)
  }
}

combineDissimEnergyPlotsInPDF = function(
  withLegend = F,
  fileIdx,
  axisFontSize = 1.5
){
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\DissimEnergyPlots\\"
  fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)
  resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"]
  resultNames = sort(resultNames)

  dissimPlotList =
    lapply(
      resultNames,
      function(rN)
        makeDissimEnergyPlot(
          results =  get(rN, envir = .GlobalEnv),
          withLegend = withLegend,
          backgroundCol = "white", textCol = "black"
        )
    )
  ggsave(
    filename = fileName,
    plot =
      grid.arrange(
        grobs = dissimPlotList,
        ncol = 1, nrow = length(resultNames),
        heights = rep(2, length(resultNames))
      ),
    height = 8.5,
    width = 4,
    dpi = 600
  )
}

makeDissimEnergyPlot = function(
  resultsList,
  printPDF = F,
  fileName,
  withLegend = F,
  textCol,
  backgroundCol,
  axisFontSize = 1.5
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
    c("#FF7F00", "#984EA3", "#377EB8", "#FFFF33", "#E41A1C", "#F781BF", "#4DAF4A")
  # c("chartreuse3", "darkorange", "blue3", "brown3", "purple", "yellow", "cyan3")
  processLabels =
    c("Nitrification", "Sulfide oxidation", "Aerobic heterotrophy", "Denitrification", "Methanogenesis", "Sulfate reduction", "Methane oxidation")
  names(processColors) = processNames
  names(processLabels) = processNames

  dissimEnergyDF$process =
    factor(
      dissimEnergyDF$process,
      levels =
        c(
          "HetMethanogenesis",
          "AutNitrif",
          "HetSulfateRed",
          "AutSulfideOxidation",
          "HetDenit",
          "HetAerobic",
          "MetMethaneOxid"
          ),
      ordered = T
    )


  processColors = processColors[names(processColors) %in% levels(dissimEnergyDF$process)]
  processLabels = processLabels[names(processLabels) %in% levels(dissimEnergyDF$process)]

  barHeights = plyr::ddply(dissimEnergyDF, "timestep", plyr::summarise, energyInJ = sum(energyInJoules))
  print(max(barHeights$energyInJ))

  dissimEnergyPlot = ggplot2::ggplot(dissimEnergyDF, ggplot2::aes(x = timestep, y = energyInJoules, fill = process)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::expand_limits(y=c(0,0.1)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = backgroundCol, colour = backgroundCol),
          panel.background = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(
      breaks = seq(1,9,1),
      labels = as.character(seq(1, 9, 1)),
      name = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 0.09, 0.03),
      labels = as.character(seq(0, 0.09, 0.03)),
      # name = "Joules"
      name = ""
    ) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = textCol, size = ggplot2::rel(axisFontSize))) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = backgroundCol),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), colour = textCol  ),
      legend.position = legendLocation) +
    ggplot2::scale_fill_manual(
      values = processColors,
      labels = processLabels
    )+
    # ggplot2::labs(x = "", y = "Joules")+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = ggplot2::rel(1.2)))


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
  axisFontSize = 1.2,
  filePrefix
){

  fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)
  resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"]
  resultNames = sort(resultNames)

  plotList =
    base::lapply(
      resultNames,
      function(rN)
        makeDissimEnergyPlot(
          results =  get(rN, envir = .GlobalEnv),
          withLegend = withLegend,
          backgroundCol = "white", textCol = "black"
        )
    )
  ggplot2::ggsave(
    filename = fileName,
    plot =
      gridExtra::grid.arrange(
        grobs = plotList,
        ncol = length(resultNames), nrow = 1,
        heights = 4
      ),
    height = 4,
    width = 8.5,
    dpi = 600
  )
}

energyBalPlot = function(resultsList, withLegend = F, backgroundCol = "white", textCol = "black", axisFontSize = 1.5){

  if(withLegend ==F) {
    legendLocation = "none"
  }else{
    legendLocation = "top"
  }

  nTS = length(resultsList)

  respEnergiesByTimeStep =
    sapply(1:nTS, function(i) sum(resultsList[[i]]$respirationEnergyVals))

  dissimEnergiesByTimeStep =
    sapply(1:nTS, function(i) sum(resultsList[[i]]$processEnergyVals$energy[resultsList[[i]]$processEnergyVals$procType == "catabolic"]))

  assimEnergiesByTimeStep =
    sapply(1:nTS, function(i) sum(resultsList[[i]]$processEnergyVals$energy[resultsList[[i]]$processEnergyVals$procType == "anabolic"]))

  energyDF =
    rbind(
      data.frame(ts = as.character(1:9), type = "Dissimilatory", energy = dissimEnergiesByTimeStep),
      data.frame(ts = as.character(1:9), type = "Assimilatory", energy = assimEnergiesByTimeStep),
      data.frame(ts = as.character(1:9), type = "Respiration", energy = respEnergiesByTimeStep)
    )

  energyDF$energyInJoules = energyDF$energy * 1000

  ggplot2::ggplot(energyDF, ggplot2::aes(x = ts, y = energyInJoules, fill = type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = c("seagreen", "darkorchid", "goldenrod")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "white")
    ) +
    ggplot2::expand_limits(y=c(-0.1,0.1)) +
    ggplot2::guides(fill = ggplot2::guide_legend(title=""))  +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = backgroundCol),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), colour = textCol  ),
      legend.position = legendLocation) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = textCol, size = ggplot2::rel(axisFontSize))) +
    # ggplot2::theme(
    #   axis.text.x = ggplot2::element_text(size = 20),
    #   axis.text.y = ggplot2::element_text(size = 20)
    # ) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0))+
    ggplot2::scale_y_continuous(
      breaks = seq(-0.1, 0.1, 0.1),
      labels = as.character(seq(-0.1, 0.1, 0.1)),
      name = ""
    )


}

#' @export
combineEnergyBalPlotsInPDF = function(
  withLegend = F,
  fileIdx,
  axisFontSize = 1.5,
  filePrefix
){

  fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)
  resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"]
  resultNames = sort(resultNames)

  nrgBalPlotList =
    base::lapply(
      resultNames,
      function(rN)
        energyBalPlot(
          results =  get(rN, envir = .GlobalEnv),
          withLegend = withLegend,
          backgroundCol = "white", textCol = "black"
        )
    )
  ggplot2::ggsave(
    filename = fileName,
    plot =
      gridExtra::grid.arrange(
        grobs = nrgBalPlotList,
        ncol = 1, nrow = length(resultNames),
        heights = rep(2, length(resultNames))
      ),
    height = 8.5,
    width = 4,
    dpi = 600
  )
}

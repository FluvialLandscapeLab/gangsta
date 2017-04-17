biomassPlot = function(resultsList, withLegend = F, backgroundCol = "white", textCol = "black"){

  if(withLegend ==F) {
    legendLocation = "none"
  }else{
    legendLocation = "top"
  }

  nTS = length(resultsList)

  hetBiomassByTimeStep =
    sapply(1:9, function(i) resultsList[[i]][[3]]$final[row.names(resultsList[[i]][[3]]) == "Het_C"])
  autBiomassByTimeStep =
    sapply(1:9, function(i) resultsList[[i]][[3]]$final[row.names(resultsList[[i]][[3]]) == "Aut_C"])
  metBiomassByTimeStep =
    sapply(1:9, function(i) resultsList[[i]][[3]]$final[row.names(resultsList[[i]][[3]]) == "Met_C"])

  bioDF =    rbind(
    data.frame(ts = as.character(1:9), type = "Het", mols = hetBiomassByTimeStep),
    data.frame(ts = as.character(1:9), type = "Aut", mols = autBiomassByTimeStep),
    data.frame(ts = as.character(1:9), type = "Met", mols = metBiomassByTimeStep)
  )
  ggplot2::ggplot(bioDF, ggplot2::aes(x = ts, y = mols, fill = type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = c("turquoise4", "tan1", "mediumvioletred")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "white")
    ) +
    ggplot2::expand_limits(y=c(0,0.16)) +
    ggplot2::guides(fill = ggplot2::guide_legend(title=""))  +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = backgroundCol),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), colour = textCol  ),
      legend.position = legendLocation) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 20),
      axis.text.y = ggplot2::element_text(size = 20)
    ) +
    ggplot2::labs(x = "", y = "") +
    # ggplot2::geom_hline(ggplot2::aes(yintercept = 0))+
    ggplot2::scale_y_continuous(
      breaks = seq(0, 0.15, 0.05),
      labels = as.character(seq(0, 0.15, 0.05)),
      name = ""
    )
}


combineBiomassPlotsInPDF = function(
  withLegend = F,
  fileIdx,
  axisFontSize = 1.5
){
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\BiomassPlots\\"
  fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)
  resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"]
  resultNames = sort(resultNames)

  biomassPlotList =
    base::lapply(
      resultNames,
      function(rN)
        biomassPlot(
          results =  get(rN, envir = .GlobalEnv),
          withLegend = withLegend,
          backgroundCol = "white", textCol = "black"
        )
    )
  ggplot2::ggsave(
    filename = fileName,
    plot =
      gridExtra::grid.arrange(
        grobs = biomassPlotList,
        ncol = 1, nrow = length(resultNames),
        heights = rep(2, length(resultNames))
      ),
    height = 8.5,
    width = 4,
    dpi = 600
  )
}

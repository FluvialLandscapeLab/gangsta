leakInPlot = function(
  resultsList,
  withLegend = F,
  backgroundCol = "white",
  textCol = "black",
  axisFontSize = 1.5,
  removeBiomass = F,
  biomassNames = c("Het", "Aut", "Met"),
  removeInfiniteCompounds = T,
  infCompoundNames = c("Ox", "Hx"),
  compoundOrder = c("O2","Ox", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "DOM", "CO2", "Het", "Aut", "Met"),
  colorVect = c(
    O2 = "#377EB8",
    Biomass = "black",
    DOM = "#E41A1C",
    Ox = "#999999",
    SO4 = "#F781BF",
    CH4 = "#4DAF4A",
    HS = "#984EA3",
    Hx = "#999999",
    NH4 = "#FF7F00",
    NO3 = "#FFFF33",
    N2 = "#A65628",
    CO2 = "#999999"
  )
){

  # decide if you want a legend
  if(withLegend ==F) {
    legendLocation = "none"
  }else{
    legendLocation = "top"
  }

  # indicate number of time steps
  nTS = length(resultsList)

  # return leak in list for each time step as a matrix
  leakIns = sapply(1:nTS, function(ts) resultsList[[ts]]$leakInCompoundVals)
  colnames(leakIns) = as.character(1:nTS)
  leakIns = as.data.frame(leakIns)

  # add a column for time zero
  leakIns$"9999" = 0
  leakIns = leakIns[, order(names(leakIns))]

  # this next bit of code is a major hack that I don't know how to get around.
  # We need a box for O2 in the reactants for the scenario wherein O2 is a
  # source.  In the particular scenario that we want to put in the paper,
  # it is a CN model with O2, Ox, and Hx as source/sinks.  Well, O2 isn't listed
  # in the gangstas list as an infinite compound because the model is a carbon
  # and nitrogen model only.  Run the following line of code if you don't
  # believe me: getGangstaAttribute(subsetGangstas(gangstasCN_Ox.O2.Hx ,
  # "class", "compound"), "InfiniteCompound").  Okay, so I am setting the O2 box
  # to be of arbitrary size in the Reactants of the gangstasCN_Ox.O2.Hx model.
  # if(identical(resultsList, resultsCN_Ox.O2.Hx)){
  #   leakIns = rbind(leakIns, O2 = rep(0.25, nTS))
  # }

  # assign column names as time step
  colnames(leakIns) = as.character(0:nTS)

  # order compound names according to compoundOrder
  inputNames = compoundOrder[compoundOrder %in% row.names(leakIns)]
  leakIns = leakIns[match(inputNames, row.names(leakIns)), ]

  # if removing infinite compounds, do so now
  if(removeInfiniteCompounds == T) {
    leakIns = leakIns[!(row.names(leakIns) %in% infCompoundNames), ]
  }

  # if removing biomass compounds, do so now
  if(removeBiomass == T) {
    leakIns = leakIns[!(row.names(leakIns) %in% biomassNames), ]
  } else{
    bioLeakIn = colSums(leakIns[(row.names(leakIns) %in% biomassNames), ])
    leakIns = leakIns[!(row.names(leakIns) %in% biomassNames), ]
    leakIns = rbind(leakIns, Biomass = bioLeakIn)
    }

  # remove compounds in color vect that are not going to be plotted and make
  # sure that the order of the color vect matches the row names in the leak in
  # matrix
  colorVect = colorVect[names(colorVect) %in% row.names(leakIns)]
  colorVect = colorVect[match(row.names(leakIns), names(colorVect))]

  # convert matrix to data frame and add column with compound name
  leakIns = as.data.frame(leakIns)
  leakIns$Compound = row.names(leakIns)
  # convert leakIns matrix to long
  leakIns = reshape2::melt(leakIns)
  colnames(leakIns) = c("Compound", "timestep", "umol")
  # make time step a character so it will plot discretely
  leakIns$timestep = as.character(leakIns$timestep)

  # make compounnd an ordered factor so that you can control the order that the
  # compounds plot within the stacked bars
  leakIns$Compound = factor(leakIns$Compound, levels = names(colorVect), ordered = T)

  # plot it!
  ggplot2::ggplot(leakIns, ggplot2::aes(x = timestep, y = umol, fill = Compound)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = colorVect) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "white")
    ) +
    ggplot2::expand_limits(y=c(0,0.6)) +
    ggplot2::guides(fill = ggplot2::guide_legend(title=""))  +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = backgroundCol),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), colour = textCol  ),
      legend.position = legendLocation) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = textCol, size = ggplot2::rel(axisFontSize)))+
    ggplot2::scale_y_continuous(
      breaks = seq(0, 0.60, 0.20),
      labels = as.character(seq(0, 0.60, 0.20))
    ) +
    # ggplot2::labs(x = "", y = expression(paste("\u03bc", "mol")))+
    ggplot2::labs(x = "", y = "")+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = ggplot2::rel(1.2)))
}



combineLeakInPlots = function(
  withLegend = F,
  fileIdx,
  axisFontSize = 1.2) {

  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\"
  fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)

  resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"]
  resultNames = sort(resultNames)

  gangstaNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,8) == "gangstas"]
  gangstaNames = sort(gangstaNames)

  plotList =
    base::lapply(
      1:length(resultNames),
      function(rNum)
        leakInPlot(
          results =  get(resultNames[rNum], envir = .GlobalEnv),
          withLegend = withLegend
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

# this next trick makes the complete leak in list plot 4 times
combineCompleteLeakInListPlots = function(
  withLegend = F,
  fileIdx,
  axisFontSize = 1.2) {

  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\"
  fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)

  resultNames = rep("resultsCONSH_Ox.Hx", 4)

  plotList =
    base::lapply(
      1:length(resultNames),
      function(rNum)
        leakInPlot(
          results =  get(resultNames[rNum], envir = .GlobalEnv),
          withLegend = withLegend
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

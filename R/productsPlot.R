productsPlot = function(
  resultsList,
  withLegend = F,
  backgroundCol = "white",
  textCol = "black",
  axisFontSize = 1.5,
  removeBiomass = F,
  biomassNames = c("Het", "Aut", "Met"),
  removeInfiniteCompounds = T,
  infCompoundNames = c("Ox", "Hx"),
  compoundOrder = c("O2", "CH4", "DOM", "SO4", "HS", "Ox",  "Hx", "NH4", "NO3", "N2", "CO2", "Het", "Aut", "Met"),
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
  ),
  fileName = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\prodPlot.pdf",
  printPDF = F


){

  # decide if you want a legend
  if(withLegend ==F) {
    legendLocation = "none"
  }else{
    legendLocation = "top"
  }

  # indicate number of time steps
  nTS = length(resultsList)

  # return products for each time step as a data frame
  products =
    do.call(rbind,
            lapply( 1:nTS, function(ts) {
              umol = resultsList[[ts]]$compoundVals[,"final"]
              Compound = row.names(resultsList[[ts]]$compoundVals)
              timestep = as.character(rep(ts, length(umol)))
              prods = data.frame(Compound, timestep, umol)
            }
            )
    )

  # convert data frame to matrix.  doing this because I'm lazy (err, efficient)
  # and want to save time and recycle code that I wrote for the leakInPlot.
  products = labdsv::matrify(products)

  # add a column for time zero
  products$"0" = 0
  products = products[, order(names(products))]


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
  #   products = rbind(products, O2 = rep(0.25, nTS))
  # }

  # assign column names as time step
  colnames(products) = as.character(0:nTS)

  # order compound names according to compoundOrder
  inputNames = compoundOrder[compoundOrder %in% row.names(products)]
  products = products[match(inputNames, row.names(products)), ]

  # if removing infinite compounds, do so now
  if(removeInfiniteCompounds == T) {
    products = products[!(row.names(products) %in% infCompoundNames), ]
  }

  # if removing biomass compounds, do so now
  # if removing biomass compounds, do so now
  if(removeBiomass == T) {
    products = products[!(row.names(products) %in% biomassNames), ]
  } else{
    bioProds = colSums(products[(row.names(products) %in% biomassNames), ])
    products = products[!(row.names(products) %in% biomassNames), ]
    products = rbind(products, Biomass = bioProds)
  }

  # remove compounds in color vect that are not going to be plotted and make
  # sure that the order of the color vect matches the row names in the leak in
  # matrix
  colorVect = colorVect[names(colorVect) %in% row.names(products)]
  colorVect = colorVect[match(row.names(products), names(colorVect))]

  # convert matrix to data frame and add column with compound name
  products = as.data.frame(products)
  products$Compound = row.names(products)
  # convert products data frame to long
  products = reshape2::melt( products )
  colnames(products) = c("Compound", "timestep", "umol")
  # make time step a character so it will plot discretely
  products$timestep = as.character(products$timestep)

  # make compounnd an ordered factor so that you can control the order that the
  # compounds plot within the stacked bars
  products$Compound = factor(products$Compound, levels = names(colorVect), ordered = T)

  barHeights = plyr::ddply(products, "timestep", plyr::summarise, umol = sum(umol))
  print(max(barHeights$umol))

  # plot it!
  prodPlot = ggplot2::ggplot(products, ggplot2::aes(x = timestep, y = umol, fill = Compound)) +
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

  if(printPDF == TRUE) {
    pdf(
      fileName,
      onefile = T,
      paper = "USr",
      # paper = "letter",
      width = 11,
      height = 4.5
    )
    print(prodPlot)
    dev.off()
  } else {
    return(prodPlot)
  }
}

combineProductsPlots = function(
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
        productsPlot(
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




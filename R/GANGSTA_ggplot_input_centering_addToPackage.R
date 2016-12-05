

gangstaSuperPlotInput = function(results = resultsCHONS_Ox.Hx, gangstas = gangstasCHONS_Ox.Hx, aggregateBio = T) {
  imperfectResultsList = simplifyDataForAlluvialPlot(results, gangstas)
  imperfectResultsList = unlist(imperfectResultsList, recursive = F)
  imperfectResultsList = lapply(
    1:length(imperfectResultsList),
    function(ts) {
      from = as.character(imperfectResultsList[[ts]]$fromPool)
      to = as.character(imperfectResultsList[[ts]]$toPool)
      element = substr(from, nchar(from), nchar(from))
      from =  substr(from, 1, nchar(from)-2)
      to = substr(to, 1, nchar(to)-2)
      mols = imperfectResultsList[[ts]]$molsTransfered
      return(data.frame(step = ts, from = from, to = to, element = element, mols = mols, stringsAsFactors = F))
    }
  )
  perfectDF = do.call(rbind, imperfectResultsList)
  if(aggregateBio == T) {
    fromIsBio = (perfectDF$from == "Het" | perfectDF$from == "Aut" | perfectDF$from == "Met")
    toIsBio = (perfectDF$to == "Het" | perfectDF$to == "Aut" | perfectDF$to == "Met")
    perfectDF$from[fromIsBio] = "Bio"
    perfectDF$to[toIsBio] = "Bio"
  }
  perfectDF = plyr::ddply(perfectDF, c("step", "from", "to", "element"), plyr::summarise, mols = sum(mols))
}


gangstaSuperPlot = function(
  perfectDF = inputDF,
  fileName = "moneyPlot_20160517.pdf",
  makePDF = TRUE,

  compoundOrder = c("Met", "Aut", "Het", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "O2"),
  # compoundOrder = c("Bio", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "O2"),
  # compoundOrder = c("O2", "CO2", "Bio", "DOM", "Hx", "Ox", "CH4", "SO4", "HS", "NH4", "NO3", "N2"),

  elementOrder = c("C", "H", "O", "N", "S"),

  # elementColor = c(C = "#FF0000", H = "#C0C0C0", O = "#000066", N = "#006400", S = "#FFFF00"),
  # elementColor = c(C = "#B22222", H = "#C0C0C0", O = "#0000CD", N = "#228B22", S = "#FFD700"),

  elementColor =
    # rep("black", 5),

    c(C = "#B22222", H = "#C0C0C0", O = "#0000CD", N = "#228B22", S = "magenta3"),


  opacity = 0.9,
  compoundSpacing = 0.2,
  compoundBoxWidth = 0.2,
  transferGap = 0.05,
  poolGap = 0.01,
  minPoolMols = 0,
  minTransferMols = 0,
  justBoxes = F,
  backgroundCol = "black",
  textCol = "white",
  axisFontSize,
  yAxisMaxMols = 0
){

  elementColor = rgb(t(col2rgb(elementColor)/255), alpha = opacity)
  names(elementColor) = elementOrder

  plotDF = perfectDF #getPlotDF(plotData)
  # plotDF = plotDF[plotDF$element %in% c("S"), ]

  # sum pools by "from" to get thickness of transfers out and by "to" to get thickness of tranfers in.
  startPoolDF = plyr::ddply(plotDF, c("step", "from", "element"), plyr::summarise, mols = sum(mols))
  endPoolDF = plyr::ddply(plotDF, c("step", "to", "element"), plyr::summarise, mols = sum(mols))
  # adjust step of transfers out; then they are transfers in for next step.
  endPoolDF$step = endPoolDF$step + 1

  # "from" and "to" are the compound names in each DF
  names(startPoolDF)[which(names(startPoolDF) %in% "from")] = "compound"
  names(endPoolDF)[which(names(endPoolDF) %in% "to")] = "compound"

  # This creates a data frame with one row for each pool for steps 1 and n, and
  # two rows for each pool (one row for incoming and one for outgoing mols) for
  # steps from 2 to n-1.
  poolDF = rbind(startPoolDF, endPoolDF)
  # The pool height should be the max of incoming or outgoing.
  poolDF$minPoolMols = minPoolMols
  poolDF$poolGap = poolGap
  poolDF = plyr::ddply(poolDF, c("step", "compound", "element"), plyr::summarise, height = max(c(minPoolMols, max(mols))))

  # add the poolGap to mols; use plusGapHeight later to calculate yVals of pools
  poolDF$plusGapHeight = poolDF$height + poolGap


  # order according to the element order list
  activeElements = elementOrder[elementOrder %in% poolDF$element]
  rowOrder = match(poolDF$element, activeElements)
  poolDF = poolDF[order(rowOrder),]
  poolDF = poolDF[order(poolDF$compound),]
  poolDF = poolDF[order(poolDF$step),]
  poolDF$priorHeight = plyr::ddply(poolDF, c("step", "compound"), plyr::summarise, priorHeight = c(0, cumsum(height))[1:length(height)])$priorHeight
  poolDF$priorHeightPlusGap = plyr::ddply(poolDF, c("step", "compound"), plyr::summarise, priorHeightPlusGap = c(0, cumsum(plusGapHeight))[1:length(plusGapHeight)])$priorHeightPlusGap

  badElements = !(unique(poolDF$element) %in% elementOrder)
  if (any(badElements)) stop("Following elements in the data are not in the element ordering vector: ", paste(unique(poolDF$element)[badElements], collapse = ", "))

  compoundDF = plyr::ddply(poolDF, c("step", "compound"), plyr::summarise, height = sum(height), plusGapHeight = sum(plusGapHeight))
  badCompounds = !(unique(compoundDF$compound) %in% compoundOrder)
  if (any(badCompounds)) stop("Following compounds in the data are not in the compound ordering vector: ", paste(unique(compoundDF$compound)[badCompounds], collapse = ", "))

  #poolDF = join(poolDF, compoundDF, by = c("step", "compound"))

  #poolDF$cmpdHeight = sapply(1:nrow(poolDF), function(.i) compoundDF[which(compoundDF$step == poolDF$step[.i] & compoundDF$compound == poolDF$compound[.i]),]$height)

  rowDF = plyr::ddply(compoundDF, c("compound"), plyr::summarise, height = max(height), plusGapHeight = max(plusGapHeight))

  #order the row height dataframe according to the inverse of the compoundOrder
  activeCompounds = compoundOrder[compoundOrder %in% rowDF$compound]
  rowOrder = length(activeCompounds) + 1 - match(rowDF$compound, activeCompounds)
  rowDF = rowDF[order(rowOrder),]
  rowDF$rowTop = cumsum(rowDF$plusGapHeight + compoundSpacing)
  rowDF$rowCenter = rowDF$rowTop - (rowDF$plusGapHeight + compoundSpacing)/2

  # begin plot rendering
  river = ggplot2::ggplot()

  # if yAxisMaxMols != 0, Adjust compound spacing to make max of y axis equal yAxisMaxMols
  if(yAxisMaxMols > 0) {
    # find current max value for y axis
    yScaleMax = max(rowDF$rowTop)
    # if current value exceeds requested value, ooops...
    if(yScaleMax > yAxisMaxMols) {
      stop("The 'yAxisMaxMols' parameter is too small.  Summed compound mols (", yScaleMax, ") can't exceed the 'yAxisMaxMols' parameter.")
    }
    # determine the height of a space that would center the whole graph on the y Axis
    centeringGap = (yAxisMaxMols - yScaleMax)/2
    rowDF$rowTop = rowDF$rowTop + centeringGap
    rowDF$rowCenter = rowDF$rowCenter + centeringGap
    river = river + ggplot2::coord_cartesian(ylim = c(0, yAxisMaxMols))
  }


  # move y values of row tops to the poolDFrowDF$plusGapHeight + compoundSpacing
  poolDF = plyr::join(poolDF, rowDF[,c("compound", "rowCenter")], by = "compound")
  compoundHeightDF = data.frame(step = compoundDF$step, compound = compoundDF$compound, compoundPlusGapHeight = compoundDF$plusGapHeight, stringsAsFactors = F)
  poolDF = plyr::join(poolDF, compoundHeightDF, by = c("step", "compound"))

  # caluelate y values for top of each pool
  poolDF$compoundTop = poolDF$rowCenter + (poolDF$compoundPlusGapHeight - poolGap)/2
  poolDF$poolTop = poolDF$compoundTop - poolDF$priorHeightPlusGap
  poolDF$y = poolDF$poolTop - poolDF$height/2

  # render pools as a ggplot
  river = river +
    ggplot2::scale_fill_manual(values = elementColor) +
    ggplot2::geom_tile(
      data = poolDF,
      mapping = ggplot2::aes(
        x=step,
        y=y,
        height = height,
        width = compoundBoxWidth,
        fill = element
        )
      )

  # sort pools by to, from, element, and step
  rowOrder = match(plotDF$to, activeCompounds)
  plotDF = plotDF[order(rowOrder),]
  rowOrder = match(plotDF$from, activeCompounds)
  plotDF = plotDF[order(rowOrder),]
  plotDF = plotDF[order(plotDF$element),]
  plotDF = plotDF[order(plotDF$step),]

  # make a cumulative sum of pool height for all rows with same step, from
  # compound, and element.  This gives the offset from the y value of the top of
  # the pool, which will be used to calcualte the y value where the transfer comes
  # in.
  plotDF$fromPrior = 0.0
  for (i in 2:nrow(plotDF)) {
    if(all(plotDF[i-1, c("step", "from", "element")] == plotDF[i, c("step", "from", "element")])) {
      plotDF[i, "fromPrior"] = plotDF[i-1, "fromPrior"] + plotDF[i-1, "mols"]
    }
  }

  # repeat two block above, but do it for outgoing transfers.
  rowOrder = match(plotDF$from, activeCompounds)
  plotDF = plotDF[order(rowOrder),]
  rowOrder = match(plotDF$to, activeCompounds)
  plotDF = plotDF[order(rowOrder),]
  plotDF = plotDF[order(plotDF$element),]
  plotDF = plotDF[order(plotDF$step),]

  plotDF$toPrior = 0.0
  for (i in 2:nrow(plotDF)) {
    if(all(plotDF[i-1, c("step", "to", "element")] == plotDF[i, c("step", "to", "element")])) {
      plotDF[i, "toPrior"] = plotDF[i-1, "toPrior"] + plotDF[i-1, "mols"]
    }
  }

  # move the pool y location data over to the plotDF, where each row represents a
  # transfer. First for "from"...
  names(poolDF)[names(poolDF) == "compound"] = "from"
  names(poolDF)[names(poolDF) == "poolTop"] = "fromTop"
  plotDF = plyr::join(plotDF, poolDF[,c("step", "from", "fromTop", "element")], by = c("step", "from", "element"))

  # then for to...
  names(poolDF)[names(poolDF) == "from"] = "to"
  names(poolDF)[names(poolDF) == "fromTop"] = "toTop"
  poolDF$step = poolDF$step - 1
  plotDF = plyr::join(plotDF, poolDF[,c("step", "to", "toTop", "element")], by = c("step", "to", "element"))

  # Calculate y values of middle of each end of transfer.
  plotDF$fromY = plotDF$fromTop - plotDF$fromPrior - plotDF$mols/2
  plotDF$toY = plotDF$toTop - plotDF$toPrior - plotDF$mols/2

  if(justBoxes == FALSE) {
    plotIt = list()
    for(.i in 1:nrow(plotDF)) {
      plotIt[[.i]] = data.frame(
        x = c((plotDF[.i, "step"] + compoundBoxWidth/2 + transferGap),(plotDF[.i, "step"] + 1 - compoundBoxWidth/2 - transferGap)),
        y = c(plotDF[.i, "fromY"], plotDF[.i, "toY"]),
        element = plotDF[.i, "element"],
        var = plotDF[.i, "mols"]/2
      )
      if(plotDF[.i, "mols"] < minTransferMols) {
        river = river + ggplot2::geom_line(data = plotIt[[.i]], mapping = ggplot2::aes(x = x, y = y, colour = element))
      } else {
        river = river + ggplot2::geom_ribbon(data = plotIt[[.i]], mapping = ggplot2::aes(
          x = x,
          # y = y,
          ymin = y - var,
          ymax = y + var,
          fill = element)
          )
      }
    }
  }

  #river = saveriver

  river = river +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = backgroundCol, colour = backgroundCol), panel.background = ggplot2::element_blank())

  y.axisDF = plyr::ddply(poolDF, "to", plyr::summarise,  y = median(y))
  y.axisDF = y.axisDF[with(y.axisDF, order(y)), ]

  river = river +
    ggplot2::scale_y_continuous(breaks = y.axisDF$y, labels = y.axisDF$to, name = "") +
    ggplot2::scale_x_continuous(breaks = seq(1.5,9.5,1), labels = as.character(seq(1, 9, 1)), name = "") +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = textCol, size = ggplot2::rel(axisFontSize))) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = backgroundCol),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.5), colour = backgroundCol),
      legend.position = "none"
      # legend.position = "top"
    )

  river = river + ggplot2::theme(plot.margin = grid::unit(c(0.5, 1, -0.5, 0.5), "lines"))
  if(makePDF == TRUE) {
    pdf(
      fileName,
      onefile = T,
      paper = "USr",
      # paper = "letter"
      width = 10.5,
      height = 7
    )
    print(river)
    dev.off()
  } else{
    return(river)
  }
  # river
}


makeFileName =
  function(
    fileID,
    filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\RiverPlots\\"
  ){
    paste0(filePrefix, fileID, ".pdf")
  }




makeOutput =
  function(
    inputList =
      list(
        CN.anoxic = list(activeElements = c("C", "N"), sourceSinks =  c("Ox", "Hx")),
        CN.oxic = list(activeElements = c("C", "N"), sourceSinks = c("Ox", "O2", "Hx")),
        CNO = list(activeElements = c("C", "O", "N"), sourceSinks = c("Ox", "Hx")),
        CHONS = list(activeElements = c("C", "O", "N", "S", "H"), sourceSinks = c("Ox", "Hx"))
      )
  ){
    lapply(inputList, function(input) CNOSH_Any(input[["activeElements"]], input[["sourceSinks"]]))
  }

makePlots =
  function(
    resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"],
    gangstaNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,8) == "gangstas"],
    pdf = F,
    aggregateBio = T,
    elementalCyclesToPlot = NULL,
    yAxisMaxMols,
    axisFontSize
  ){
    resultNames = sort(resultNames)
    gangstaNames = sort(gangstaNames)

    fileIDXs = resultNames

    gangstaPerfectDFList =
      lapply(1:length(resultNames), function(i)
        perfectDF = gangstaSuperPlotInput(
          results = get(resultNames[i], envir = .GlobalEnv),
          gangstas = get(gangstaNames[i], envir = .GlobalEnv),
          aggregateBio = aggregateBio
        )
      )

    if(!is.null(elementalCyclesToPlot)) {
      cycleNames = do.call(paste0, as.list(elementalCyclesToPlot))
      fileIDXs  = paste0(resultNames, "_", cycleNames)
      gangstaPerfectDFList =
        lapply(
          gangstaPerfectDFList, function(DF)
            DF = DF[DF$element %in% elementalCyclesToPlot, ]
        )
    }

    plotList =
      lapply(
        1:length(resultNames), function(i)
          print(
            gangstaSuperPlot(
              perfectDF = gangstaPerfectDFList[[i]],
              makePDF = pdf,
              fileName = makeFileName(fileID = fileIDXs[i]),
              backgroundCol = "white", textCol = "black",
              axisFontSize = axisFontSize,
              yAxisMaxMols = yAxisMaxMols
            )
          )
      )
    if(pdf == F) return(plotList)
  }


combineRiverPlotsInPDF = function(
  withLegend = F,
  fileIdx,
  resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"],
  gangstaNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,8) == "gangstas"],
  aggregateBio = F,
  elementalCyclesToPlot = NULL,
  cyclesSeparate = T,
  axisFontSize,
  yAxisMaxMols = 0
){
  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\RiverPlots\\"
  fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)

  resultNames = sort(resultNames)
  gangstaNames = sort(gangstaNames)

  numColumns =
    ifelse(
      (length(elementalCyclesToPlot) > 1 && cyclesSeparate),
      length(elementalCyclesToPlot),
      1
    )
  numRows = length(resultNames)

  gangstaPerfectDFList =
    lapply(1:length(resultNames), function(i)
      perfectDF = gangstaSuperPlotInput(
        results = get(resultNames[i], envir = .GlobalEnv),
        gangstas = get(gangstaNames[i], envir = .GlobalEnv),
        aggregateBio = aggregateBio
      )
    )

  if(length(elementalCyclesToPlot)>1) {
    if(cyclesSeparate == F){
      gangstaPerfectDFList =
        lapply(gangstaPerfectDFList, function(DF) DF[DF$element %in% elementalCyclesToPlot, ])
    } else {
      numOfElCycles = length(elementalCyclesToPlot)
      gangstaPerfectDFListOrig = gangstaPerfectDFList
      gangstaPerfectDFList = list()
      for(i in 1:length(gangstaPerfectDFListOrig)) {
        tempDF = gangstaPerfectDFListOrig[[i]]
        tempDFList = lapply(elementalCyclesToPlot, function(el) tempDF[tempDF$element == el, ])
        gangstaPerfectDFList = append(gangstaPerfectDFList, tempDFList, length(gangstaPerfectDFList))
      }
    }
  }
  # https://github.com/wch/ggplot2/wiki/New-theme-system
  new_theme_empty <- ggplot2::theme_bw()
  new_theme_empty$line <- ggplot2::element_blank()
  new_theme_empty$rect <- ggplot2::element_blank()
  new_theme_empty$strip.text <- ggplot2::element_blank()
  new_theme_empty$axis.text <- ggplot2::element_blank()
  new_theme_empty$plot.title <- ggplot2::element_blank()
  new_theme_empty$axis.title <- ggplot2::element_blank()
  new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")
  ggplot2::ggsave(filename = fileName,
                  gridExtra::grid.arrange(
                    grobs =
                      lapply(
                        gangstaPerfectDFList,
                        function(DF)
                          if(nrow(DF) > 0) {
                            gangstaSuperPlot(
                              DF,
                              makePDF = F,
                              backgroundCol = "white",
                              textCol = "black",
                              axisFontSize = axisFontSize,
                              yAxisMaxMols = yAxisMaxMols
                            )
                          } else {
                            ggplot2::ggplot(mapping = ggplot2::aes(a, b), data = data.frame(a = 0, b = 0)) + new_theme_empty
                          }
                      ),
                    ncol = numColumns, nrow = numRows
                    # ,
                    # heights = rep(2, length(riverPlotList))
                  ),
                  width = 11,
                  height = 8.5,
                  units = "in",
                  dpi = 600
  )
}


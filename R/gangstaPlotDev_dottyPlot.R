plotIt = function(results, dotFactor, lineFactor, tag, excludeCompounds = character(0)) {

  # get the compound values
  compoundList = lapply(results, "[[", "compoundVals")
  compoundList = lapply(compoundList, function(x) subset(x, !(row.names(x) %in% excludeCompounds)))
  # get the pool vals
  poolDataList = lapply(results, "[[", "poolVals")
  # get transfer vals
  transferList = lapply(results, "[[", "molTransferVals")
  # get energy vals
  energyValList = lapply(results, "[[", "processEnergyVals")
  energyValList = lapply(energyValList, function(x) subset(x, x$procType == "catabolic"))
  newRowNames = c(
    HetAerobic = "Aerobic Heterotrophy",
    HetDenit = "Denitrification",
    HetSulfateRed = "Sulfate Reduction",
    HetMethanogenesis = "Methanogensis",
    AutNitrif = "Nitrification",
    AutSulfideOxidation = "Sulfide Oxidation",
    MetMethaneOxid = "Methane Oxidation")
  energyValList =
    lapply(
      energyValList,
      function(x) {
        row.names(x) = newRowNames[row.names(x)]
        return(x)
      }
    )

  dotDataList = list(compoundList, energyValList, poolDataList)
  transferDataList = list(list(), list(), transferList)
  variables = list(c("initial", "final"), c("energy"), c("initial", "final"))
  dotColors = list(c("red", "black"), c("black"), c("red", "black"))
  graphTitles = list(
    expression(paste("(a) Available reagents (", mu, "mols)")),
    "(b) Metabolic energy yield (kJ)",
    expression(paste("(c) Elemental pools and transfers (", mu, "mols)"))
  )
  xOffset = list(0,0.5,0)

  layout(matrix(c(1,2,3), 3, 1), heights = sapply(dotDataList, function(x) nrow(x[[1]])) + 1)
  par(oma = c(6,11,3,1))
  mapply(subplotIt, dotDataList, transferDataList, variables, dotColors, xOffset, graphTitles, dotFactor, lineFactor, tag)
  axis(side = 1, at = 0:length(transferList), labels = F, cex.axis = 1.0)
  axis(side = 1, at = 0:(length(transferList)-1) + 0.5, labels = 1:length(transferList), cex.axis = 1.0, tick = F)
  mtext("Time Step", 1, 3, cex = 0.8)
  titleElements = paste(unlist(strsplit(tag, split ="")), collapse = ", ")
  titleSourceSinks = paste(excludeCompounds, collapse = ", ")
  titleTag = paste("Elements:", titleElements, "; Source/sinks:", titleSourceSinks)
  title(main=titleTag, outer=TRUE, cex.main = 1.5)
}

subplotIt = function(dotList, transferList, variables, dotColors, xOffset, graphTitle, dotFactor, lineFactor, tag) {

  # make Y values for each  and add  names
  yVals = nrow(dotList[[1]]):1
  names(yVals) = row.names(dotList[[1]])

  # add the iteration and rownumdber as a column to each Val dataframe
  dotList = mapply(function(itr, pv) cbind(itr, yVals, pv), 1:length(dotList), dotList, SIMPLIFY = F)
  # rbind into a single data frame
  dotValData = do.call(rbind, dotList)

  if(length(transferList)>0) {
    # add iteration number to each transferVal dataframe
    transferList =  mapply(function(itr, tv) cbind(itr, tv), 1:length(dotList), transferList, SIMPLIFY = F)
    # extract transfers that are non-zero
    transferList = lapply(transferList, function(x) subset(x, subset = (x$molTransfered > 0.0)))
    # rbind into a single data frame
    transferValData = do.call(rbind, transferList)
    lineXvalList = lapply(transferValData$itr, function(x) c(x-1, x))
    lineYvalList = mapply(function(f, t) c(yVals[f], yVals[t]), as.character(transferValData$from), as.character(transferValData$to), SIMPLIFY = F)
  }

  maxItr = max(dotValData$itr)
  maxYVal = max(yVals)
  extraY = 1/sqrt(maxYVal)

  par(mar = c(0.8, 0.5, 1.5, 0.5))
  plot(1, type="n",
       xlim = c(0, maxItr),
       ylim = c(1-extraY, maxYVal+extraY),
#       ylab = "",
       xaxt = "n",
       yaxt = "n"
  )
  par(yaxp = c(1, maxYVal, max(1,maxYVal-1)))
  grid(nx = 0, ny = NULL)
  if(length(transferList)>0) {
    mapply(lines, x = lineXvalList, y = lineYvalList, lwd = transferValData$molTransfered * lineFactor, col = "darkgrey")
  }

  for(i in 1:length(variables) ) {
    points(dotValData$itr+i-2+xOffset,
           dotValData$yVals,
           pch = 16,
           cex = radiusOfArea(dotValData[,variables[i]], dotFactor),
           col = dotColors[i]
    )
  }

  axis(side = 2, at = yVals, labels = names(yVals), las = 2, cex.axis = 1.0)
  axis(side = 1, at = 0:maxItr, labels = F)

  mtext(graphTitle, at = -0.3, line = 0.2, adj = 0)
#  xaxp = c(0, maxItr, maxItr),

}

radiusOfArea = function(areas, radiusFactor = 1) {
    return(radiusFactor * sqrt(areas / pi))
}


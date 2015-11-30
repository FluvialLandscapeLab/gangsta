###  Imports gangsta lp file into R
readGangsta.lp = function(lpFile = file.choose()){
  lp.bgc = lpSolveAPI::read.lp(lpFile, type = "lp", verbose = "full")
  return(lp.bgc)
}

### This function allows the user to access the columns of the lp object.
getModelCol.lp = function(lpObject, colnum) {
  sapply(1:dim(lpObject)[1], function(i) get.mat(lpObject, i, colnum))
}

###  This function allows the user to view the values of the lp object.
showMatrixVals.lp = function(varName, matrix = lpMatrix, lpObject) {
  varCol = which(colnames(matrix) == varName)
  varRows = which(getModelCol(lpObject, varCol) != 0)
  matrix[varRows, varCol]
}

###  This block of code builds a matrix of decision variables and constraints
viewMatrix.lp = function(lpObject){
  constraintVarNames = dimnames(lpObject)[[1]]
  decisionVarNames = dimnames(lpObject)[[2]]
  nvars = length(decisionVarNames)
  nconstraints = length(constraintVarNames)
  lpMatrix = matrix(mapply(function(i,j) lpSolveAPI::get.mat(lpObject, i, j),
                           i=1:nconstraints, j=rep(1:nvars, each=nconstraints)),
                    nrow = nconstraints)
  colnames(lpMatrix) = decisionVarNames
  rownames(lpMatrix) = constraintVarNames
  return(lpMatrix)
}

### View output in data frame
solvedDataFrame.lp = function(lpObject, simple = TRUE) {
  variableOutput = lpSolveAPI::get.variables(lpObject)[order(dimnames(lpObject)[[2]])]
  names(variableOutput) = dimnames(lpObject)[[2]][order(dimnames(lpObject)[[2]])]
  variableOutput = data.frame(variableOutput)
  if(simple == FALSE){
    return(variableOutput)
  } else {
    variableOutputSimple = subset(variableOutput, variableOutput != 0)
    variableOutputSimple[order(row.names(variableOutputSimple)),]
    return(variableOutputSimple)
  }
}

compoundDifs = function(gangstaObjects, lpObject, simple = F) {
  compounds = subsetGangstas(gangstaObjects, "class", gangstaClassName("comp"))
  compoundNames = getGangstaAttribute(compounds, gangstaAttributeName("name"))

  initCompoundNames = makeCompoundStartMassVars(compoundNames)
  finalCompoundNames = makeCompoundEndMassVars(compoundNames)

  df.lp = solvedDataFrame.lp(lpObject, simple = F)

  startVals = df.lp[initCompoundNames,]
  endVals = df.lp[finalCompoundNames,]
  allDF = data.frame(change = endVals - startVals, initial = startVals, final = endVals)
  row.names(allDF) = compoundNames

  if(simple) allDF = subset(allDF, (allDF$initial != 0) | (allDF$final != 0))
  return(allDF)
}


### Make a data frame with three columns: start values, end values, and changes in values
poolDifs = function(gangstaObjects, lpObject, simple = T) {
  pools = subsetGangstas(gangstaObjects, "class", gangstaClassName("pool"))
  poolNames = getGangstaAttribute(pools, gangstaAttributeName("name"))
  elementName = getGangstaAttribute(pools, gangstaAttributeName("element"))
  elementIdx = order(elementName)

  initPoolNames = makePoolStartMassVars(poolNames)
  finalPoolNames = makePoolEndMassVars(poolNames)

  df.lp = solvedDataFrame.lp(lpObject, simple = F)

  startVals = df.lp[initPoolNames,]
  endVals = df.lp[finalPoolNames,]
  allDF = data.frame(change = endVals - startVals, initial = startVals, final = endVals)
  row.names(allDF) = poolNames

  allDF = allDF[elementIdx, ]

  if(simple) allDF = subset(allDF, (allDF$initial != 0) | (allDF$final != 0))
  #  col.names(appDF) = c("initial", "final", "change")
  return(allDF)
}

### Make a vector of mass transfers between pools
massTransfers = function(gangstaObjects, lpObject, simple = FALSE, byProcess = FALSE) {

  ### get classes, names, attributes, and gangstaObjectss that are used later in this function
  transClassName = gangstaClassName("trans")
  nameAttrName = gangstaAttributeName("name")

  transformations = subsetGangstas(gangstaObjects, "class", transClassName)

  ### make metabolic transformation variable names from gangstaObjects
  transformationNames = getGangstaAttribute(transformations, nameAttrName)
  transformationMassTransVars = makeTransformationMassTransVars(transformationNames)

  ### make two named vectors: one of the from pools and one of the to pools
  massTransferFroms = unlist(lapply(transformations, "[[", "from") )
  massTransferTos = unlist(lapply(transformations, "[[", "to") )
  names(massTransferFroms) = transformationMassTransVars
  names(massTransferTos) = transformationMassTransVars

  ### make a data frame of lpsolve results that will be subsetted using the transformationMassTransVars from gangsta
  massTransferDF = solvedDataFrame.lp(lpObject, simple = FALSE)
  subsettingCriteria = row.names(massTransferDF) %in% transformationMassTransVars
  massTransferNames = row.names(massTransferDF)[subsettingCriteria]
  if(length(transformationMassTransVars) != length(massTransferNames))
    stop("The number of mass transfer variables generated from your gangsta list doesn't match the number of mass transfer variable names in your lpSolve file.")
  if(any(transformationMassTransVars %in% massTransferNames) == FALSE)
    stop("Your names of mass transfers in your gangsta list doesn't match your mass transfer variable names in your lpSolve file.")

  ### ensure that the order of the froms and tos matches that of the transfers in the dataframe
  massTransferFroms = massTransferFroms[match(massTransferNames, transformationMassTransVars)]
  massTransferTos = massTransferTos[match(massTransferNames, transformationMassTransVars)]

  ### subset the data frame
  massTransfers = massTransferDF[subsettingCriteria,]
  massTransfers = data.frame(massTransferFroms, massTransferTos, massTransfers)
  names(massTransfers) = c("fromPool", "toPool", "massTransfered")
  row.names(massTransfers) = massTransferNames

  if(!byProcess) {
    massTransfers = plyr::ddply(massTransfers, c("fromPool", "toPool"), plyr::summarise,
                                massTransfered = sum(massTransfered))
  }

  if(simple) {
    massTransfers = subset(massTransfers, (massTransfers$massTransfered != 0))
  }

  return(massTransfers)
}

### get energy produced/consumed by processes

processEnergies = function(gangstaObjects, lpObject,
                           condense = TRUE,
                           simple = FALSE,
                           catabolic = TRUE,
                           anabolic = TRUE,
                           decayString = "Decay",
                           stringsToStripFromEnergyVarNames = c("ofDOM", "ofBiomass")) {

  if(!any(c(catabolic, anabolic))) {
    stop("You can request anabolic, catabolic, or all energy variables.  Request for neither has failed.")
  }

  processes = subsetGangstas(gangstaObjects, "class", gangstaClassName("proc"))
  processNames = getGangstaAttribute(processes, gangstaAttributeName("name"))
  processEnergy = getGangstaAttribute(processes, gangstaAttributeName("energy"))

  ## Decay is neither catabolic or anabolic because it neither generates nor consumes energy, so we remove it:
  processEnergy = processEnergy[-grep(decayString, names(processEnergy))]
  processNames = processNames[-grep(decayString, processNames)]

  simplifiedProcessNames = processNames

  for(stringToStrip in stringsToStripFromEnergyVarNames) {
    stringsToRemoveId = grep(stringToStrip, simplifiedProcessNames)
    simplifiedProcessNames[stringsToRemoveId] =
      sub(stringToStrip, "", simplifiedProcessNames)[stringsToRemoveId]
  }
  names(simplifiedProcessNames) = processNames

  catabolicNames = names(processEnergy)[processEnergy>0]
  anabolicNames = names(processEnergy)[processEnergy<0]

  procEnergyVarNames = makeProcessEnergyVars(processNames)

  df.lp = solvedDataFrame.lp(lpObject, simple = F)

  energyVals = df.lp[procEnergyVarNames,]
  names(energyVals) = processNames

  energyVals = sort(energyVals)
  simplifiedProcessNames = simplifiedProcessNames[match(names(energyVals), names(simplifiedProcessNames))]

  anabolicVals = energyVals[match(anabolicNames, names(energyVals), nomatch = 0)]
  catabolicVals = energyVals[match(catabolicNames, names(energyVals), nomatch = 0)]

  anabolicSimpleNames = simplifiedProcessNames[match(anabolicNames, names(energyVals), nomatch = 0)]
  catabolicSimpleNames = simplifiedProcessNames[match(catabolicNames, names(energyVals), nomatch = 0)]

  anabolicDF = data.frame(energy = anabolicVals, procSimple = anabolicSimpleNames,  row.names = names(anabolicVals), procType = rep("anabolic", length(anabolicSimpleNames)))
  catabolicDF = data.frame(energy = catabolicVals, procSimple = catabolicSimpleNames, row.names = names(catabolicVals), procType = rep("catabolic", length(catabolicSimpleNames)))

  if(catabolic) {
    returnDF = catabolicDF
    if(anabolic) {
      returnDF = rbind(catabolicDF, anabolicDF)
    }
  } else {
    returnDF = anabolicDF
  }
  if(condense) {
    returnDF$returnDFSorter = c(1, rep(0, nrow(returnDF)-1))
    for(i in 2:nrow(returnDF)) {
      returnDF$returnDFSorter[i] = ifelse((returnDF$procSimple[i-1] == returnDF$procSimple[i]),
                                          returnDF$returnDFSorter[i-1],
                                          returnDF$returnDFSorter[i-1]+1)
    }
    returnDF = plyr::ddply(.data = returnDF,
                           .variables = c("returnDFSorter", "procSimple", "procType"),
                           .fun = "summarise",
                           energy = sum(energy)
    )
    names(returnDF)[match("summarise", names(returnDF),  nomatch = 0)] = "energy"
    returnDF = returnDF[, c("procSimple", "energy", "procType")]
  }
  if(simple) returnDF = subset(returnDF, energy != 0)
  return(returnDF)
}

### Make plots of multiple lp models
#############  THIS NEEDS UPDATING WHEN YOU FIGURE OUT HOW TO MAKE SEQUENTIAL UPDATES TO COMPOUNDS
massTransfersPlot = function(gangstaObjects,
                             lpResultsList,

                             tag = tag,
                             sourceSinks = sourceSinks,

                             lineMult = 24,
                             dotMult = 6,
                             energyPtMultiplier = 1000,

                             ###infrequently changed arguments:
                             xBound = 0.25,
                             nPlots = ((2*length(lpResultsList)) - 1),
                             layoutMatrix = matrix(nrow = 1, ncol = nPlots, data = 1:nPlots),
                             layoutWidths = c(2, rep(c(1,2), length(lpResultsList)-1)),
                             layoutMarginsFirstPlot = c(1,0,1,0),
                             layoutMarginsSubsequentPlot = c(1,0,1,0),
                             layoutMarginsLastPlot = c(1,0,1,0),
                             parOuterMarginSettings = c(1.75,11,2,11)
){
  numberOfIterations = length(lpResultsList)


  titleElements = paste(unlist(strsplit(tag, split ="")), collapse = ", ")
  titleSourceSinks = paste(sourceSinks, collapse = ", ")
  titleTag = paste("Elements:", titleElements, "; Source/sinks:", titleSourceSinks)

  ### need an error check that goes here to ensure that the names of each of the results
  ### items has the same sets of row names for the pool diffs data frame and the same from/to columns
  ### for the mass transfers data frame
  numOfPools = nrow(lpResultsList[[1]]$poolVals)
  numOfTransfers = nrow(lpResultsList[[1]]$massTransferVals)

  ### Determine pool locations in results list and determine where y values are plotted
  poolNames =
    row.names(lpResultsList[[1]]$poolVals)
  poolYVal =
    seq(numOfPools, 1, -1)
  names(poolYVal) = poolNames
  massTransferStartPoolIndex =
    match(lpResultsList[[1]]$massTransferVals$fromPool, poolNames)
  massTransferEndPoolIndex =
    match(lpResultsList[[1]]$massTransferVals$toPool, poolNames)
  arrowStartYVal = poolYVal[massTransferStartPoolIndex]
  arrowEndYVal = poolYVal[massTransferEndPoolIndex]

  ###  Specify x and y locations of points in mass transfer plots
  xLocs = rep(c(-xBound/2, xBound/2), each = numOfPools)
  yLocs = rep(poolYVal, 2)

  # dev.new(height = 5, width = 9)
  layout(layoutMatrix, layoutWidths)
  par(oma = parOuterMarginSettings)

  ### Set up margins of first, middle, and last plots
  marginList = list(layoutMarginsFirstPlot)
  if(numberOfIterations>1) {
    marginList = c(marginList, rep(list(layoutMarginsSubsequentPlot), numberOfIterations-2), list(layoutMarginsLastPlot))
  }

  #####################  Begin plotting  #####################
  for(i in 1:numberOfIterations){
    ### Set margins, build plot
    par(mar = marginList[[i]])

    ### Calculate point sizes in mass transfer plots
    initialVal = ifelse(lpResultsList[[i]]$poolVals$initial>0, lpResultsList[[i]]$poolVals$initial, 0)
    finalVal = ifelse(lpResultsList[[i]]$poolVals$final>0, lpResultsList[[i]]$poolVals$final, 0)
    pointSizes =
      dotMult * sqrt(c(initialVal, finalVal) / pi)
    plot(yLocs ~ xLocs,
         type = "n",
         yaxt = "n", ylab = "",
         xaxt = "n", xlab = "", xlim = c(-xBound, xBound), xaxs = "i"
    )
    points(xLocs, yLocs, pch = 16, cex = pointSizes)
    axis(1,0,labels = i, cex.axis = 1.5)
    ### give first and last plots y-axes
    if (i == 1) axis(side = 2, at = poolYVal, labels = poolNames, las = 2, cex.axis = 1.5)
    if (i == numberOfIterations) axis(side = 4, at = poolYVal, labels = poolNames, las = 2, cex.axis = 1.5)

    ### Calculate the line weights for mass transfers. dependent on par()$cex so
    ### must be calculated immediately before plot.
    lineWeights = lineMult * par()$cex * sqrt(lpResultsList[[i]]$massTransferVals$massTransfered / pi)
    ### add lines that show transfers between pools
    activeTransferIndex = which(lineWeights > 0)
    arrows(rep(-xBound/2, length(activeTransferIndex)),
           arrowStartYVal[activeTransferIndex],
           rep(xBound/2, length(activeTransferIndex)),
           arrowEndYVal[activeTransferIndex],
           lwd = lineWeights[activeTransferIndex],
           length = 0)

    ### Make leak in plot AND initialize values for next plot
    if(i < numberOfIterations) {
      par(mar = layoutMarginsSubsequentPlot)
      leakYLocs = poolYVal[names(lpResultsList[[i+1]]$leakInPoolVals)]
      pointSizes =
        dotMult * sqrt(lpResultsList[[i+1]]$leakInPoolVals / pi)

      plot(poolYVal ~ rep(0, length(poolYVal)),
           type = "n",
           yaxt = "n", ylab = "",
           xaxt = "n", xlab = "", xlim = c(-xBound/2, xBound/2), xaxs = "i", bty = "n"
      )
      points(rep(0,length(leakYLocs)), leakYLocs, pch = 16, cex = pointSizes)
    }
  }
  title(main=titleTag, outer=TRUE, cex.main = 1.5)


  ####  Make plots of catabolic energy of each timestep
  catabolicSubset = lpResultsList[[1]]$processEnergyVals$procType == "catabolic"
  energyDF = lpResultsList[[1]]$processEnergyVals
  numOfCatabProc = nrow(energyDF[catabolicSubset,])
  energyXLocs = rep(0, numOfCatabProc)
  energyYLocs = seq(numOfCatabProc, 1, -1)
  procNames = energyDF$procSimple[catabolicSubset]

  ### This next bit is a brittle hack-fest that I wrote for the paper.
  ### I'm happy to clean this up later if we want to add a fancy name
  ### to the gangsta objects that we can call on when we go to plot...
  procNames = ifelse(procNames == "HetAerobic", "Aerobic resp.",
                     ifelse(procNames =="HetDenit", "Denitrification",
                            ifelse(procNames == "HetSulfateRed", "Sulfate reduct.",
                                   ifelse(procNames == "HetMethanogenesis", "Methanogenesis",
                                          ifelse(procNames == "AutNitrif", "Nitrification",
                                                 ifelse(procNames == "MetMethaneOxid", "Methane oxid.",
                                                        "Oops. I should write less brittle code."))))))
  # dev.new(height = 5, width = 9)
  layout(layoutMatrix, layoutWidths)
  par(oma = parOuterMarginSettings)

  for(i in 1:numberOfIterations){

    energyDF =
      lpResultsList[[i]]$processEnergyVals[lpResultsList[[i]]$processEnergyVals$procType == "catabolic",]
    energyPtSizes =
      energyPtMultiplier * sqrt(energyDF$energy / pi)

    ### Set margins, build plot
    par(mar = marginList[[i]])
    plot(energyYLocs ~ energyXLocs,
         type = "n",
         yaxt = "n", ylab = "",
         xaxt = "n", xlab = "", xlim = c(-xBound, xBound), ylim = c(min(energyYLocs)-0.25, max(energyYLocs)+0.25), xaxs = "i"
    )
    points(energyXLocs, energyYLocs, pch = 16, cex = energyPtSizes)
    axis(1,0,labels = i, cex.axis = 1.5)
    ### give first and last plots y-axes
    if (i == 1) axis(side = 2, at = energyYLocs, labels = procNames, las = 2, cex.axis = 1.5)
    if (i == numberOfIterations) axis(side = 4, at = energyYLocs, labels = procNames, las = 2, cex.axis = 1.5)

    ### make a blank plot in between the energypoint plots
    if(i < numberOfIterations) {
      par(mar = layoutMarginsSubsequentPlot)
      plot(energyYLocs ~ energyXLocs,
           type = "n",
           yaxt = "n", ylab = "",
           xaxt = "n", xlab = "", xlim = c(-xBound/2, xBound/2), xaxs = "i", bty = "n"
      )
    }
  }
  title(main=titleTag, outer=TRUE, cex.main = 1.5)


  return("I am awesome-o.")
}

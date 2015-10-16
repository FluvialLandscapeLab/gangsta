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

### Make a data frame with three columns: start values, end values, and changes in values
poolDifs = function(gangstaObjects, lpObject, simple = T) {
  pools = subsetGangstas(gangstaObjects, "class", gangstaClassName("pool"))
  poolNames = getGangstaAttribute(pools, gangstaAttributeName("name"))
  elementName = getGangstaAttribute(pools, gangstaAttributeName("element"))
  elementIdx = order(elementName)

  initPoolNames = paste0(poolNames, ".", gangstaVarName("startSuffix"))
  finalPoolNames = paste0(poolNames, ".", gangstaVarName("endSuffix"))

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
  transClassName =      gangstaClassName("trans")
  nameAttrName =      gangstaAttributeName("name")

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
  if(byProcess == TRUE) {
    return(massTransfers)
  } else {
    massTransfers = plyr::ddply(massTransfers, c("fromPool", "toPool"), plyr::summarise,
                                massTransfered = sum(massTransfered)
    )
    return(massTransfers)
  }
}

### Make plots of multiple lp models
#############  THIS NEEDS UPDATING WHEN YOU FIGURE OUT HOW TO MAKE SEQUENTIAL UPDATES TO COMPOUNDS
massTransfersPlot = function(gangstaObjects,
                             lpResultsList = gangstaResults,
                             lineMult = 4,
                             dotMult = 10,
                             ###unfrequently changed arguments:
                             xBound = 0.25,
                             #                             xPositions = c(-0.125, 0.125),
                             nPlots = ((2*length(lpResultsList)) - 1),
                             layoutMatrix = matrix(nrow = 1, ncol = nPlots, data = 1:nPlots),
                             layoutWidths = c(2, rep(c(1,2), length(lpResultsList)-1)),
                             layoutMarginsFirstPlot = c(1,1,1,0), #c(3,12,1,0),
                             layoutMarginsSubsequentPlot = c(1,0,1,0), #c(3,0,1,0.08),
                             layoutMarginsLastPlot = c(1,0,1,1),#c(3,0,1,12),
                             parOuterMarginSettings = c(1,5,1,5)
){
  numberOfIterations = length(lpResultsList)

  layout(layoutMatrix, layoutWidths)
  par(oma = parOuterMarginSettings)

  ### need an error check that goes here to ensure that the names of each of the results
  ### items has the same sets of row names for the pool diffs data frame and the same from/to columns
  ### for the mass transfers data frame
  numOfPools = nrow(lpResultsList[[1]]$poolVals)
  numOfTransfers = nrow(lpResultsList[[1]]$massTransferVals)

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

  xLocs = rep(c(-xBound/2, xBound/2), each = numOfPools)
  yLocs = rep(poolYVal, 2)

  marginList = list(layoutMarginsFirstPlot)
  if(numberOfIterations>1) {
    marginList = c(marginList, rep(list(layoutMarginsSubsequentPlot), numberOfIterations-2), list(layoutMarginsLastPlot))
  }

  for(i in 1:numberOfIterations){
    pointSizes =
      sqrt(dotMult * c(lpResultsList[[i]]$poolVals$initial, lpResultsList[[i]]$poolVals$final) / pi)

    lineWeights = lineMult * lpResultsList[[i]]$massTransferVals$massTransfered
    par(mar = marginList[[i]])
    plot(yLocs ~ xLocs,
         type = "n",
         yaxt = "n", ylab = "",
         xaxt = "n", xlab = "", xlim = c(-xBound, xBound), xaxs = "i"
    )
    points(xLocs, yLocs, pch = 16, cex = pointSizes)
    axis(1,0,labels = i, cex.axis = 1.5)
    if (i == 1) axis(side = 2, at = poolYVal, labels = poolNames, las = 2, cex.axis = 1.5)
    if (i == numberOfIterations) axis(side = 4, at = poolYVal, labels = poolNames, las = 2, cex.axis = 1.5)

    activeTransfers = which(lineWeights > 0)
    arrows(rep(-xBound/2, length(activeTransfers)),
           arrowStartYVal[activeTransfers],
           rep(xBound/2, length(activeTransfers)),
           arrowEndYVal[activeTransfers],
           lwd = lineWeights, length = 0)
    if(i < numberOfIterations) {
      par(mar = layoutMarginsSubsequentPlot)
      leakYLocs = poolYVal[names(lpResultsList[[i+1]]$leakInVals)]
      pointSizes =
        sqrt(dotMult * lpResultsList[[i+1]]$leakInVals / pi)

      plot(poolYVal ~ rep(0, length(poolYVal)),
           type = "n",
           yaxt = "n", ylab = "",
           xaxt = "n", xlab = "", xlim = c(-xBound/2, xBound/2), xaxs = "i"
      )
      points(rep(0,length(leakYLocs)), leakYLocs, pch = 16, cex = pointSizes)
    }
  }
  return("I am awesome-o.")
}

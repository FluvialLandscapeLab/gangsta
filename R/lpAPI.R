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
  metabolicClassName =  gangstaClassName("metab")
  transClassName =      gangstaClassName("trans")

  nameAttrName =      gangstaAttributeName("name")
  procNameAttrName =  gangstaAttributeName("procName")

  transSuffix = gangstaVarName("transSuffix")

  metabolics =      subsetGangstas(gangstaObjects, "class", metabolicClassName)
  transformations = subsetGangstas(gangstaObjects, "class", transClassName)

  metabolicNames = getGangstaAttribute(metabolics, nameAttrName)

  ### make metabolic transformation variable names from gangstaObjects
  metabolicTransformations = lapply(metabolicNames, subsetGangstas, gangstaObjects = transformations, attributeName = procNameAttrName)
  metabolicTransformations = unlist(metabolicTransformations, recursive = FALSE)
  makeTransformationMassTransVars = function(transformationNames) {
    return(paste0(transformationNames, ".", transSuffix))
  }
  metabolicTransformationNames = getGangstaAttribute(metabolicTransformations, nameAttrName)
  metabolicTransformationMassTransVars = makeTransformationMassTransVars(metabolicTransformationNames)

  ### make two named vectors: one of the from pools and one of the to pools
  massTransferFroms = unlist(lapply(metabolicTransformations, "[[", "from") )
  massTransferTos = unlist(lapply(metabolicTransformations, "[[", "to") )
  massTransferFromsTosNames = paste(unlist(lapply(metabolicTransformations, "[[", nameAttrName)), transSuffix, sep = ".")
  names(massTransferFroms) = names(massTransferTos) = massTransferFromsTosNames

  ### make a data frame of lpsolve results that will be subsetted using the metabolicTransformationMassTransVars from gangsta
  massTransferDF = solvedDataFrame.lp(lpObject, simple = FALSE)
  subsettingCriteria = row.names(massTransferDF) %in% metabolicTransformationMassTransVars
  massTransferNames = row.names(massTransferDF)[subsettingCriteria]
  if(length(massTransferFromsTosNames) != length(massTransferNames))
    stop("The number of mass transfer variables generated from your gangsta list doesn't match the number of mass transfer variable names in your lpSolve file.")
  if(any(massTransferFromsTosNames %in% massTransferNames) == FALSE)
    stop("Your names of mass transfers in your gangsta list doesn't match your mass transfer variable names in your lpSolve file.")

  ### ensure that the order of the froms and tos matches that of the transfers in the dataframe
  massTransferFroms = massTransferFroms[match(massTransferNames, massTransferFromsTosNames)]
  massTransferTos = massTransferTos[match(massTransferNames, massTransferFromsTosNames)]

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
massTransfersPlot = function(gangstaObjects, lpObjects, lineMult = 1, dotMult = 1){

  nPlots = (2*length(lpObjects)) - 1
  layoutMatrix = matrix(nrow = 1, ncol = nPlots, data = 1:nPlots)
  layoutWidths = c(2, rep(c(1,2), length(lpObjects)-1))
  layout(layoutMatrix, layoutWidths)
  layoutMarginsFirstPlot = c(3,12,1,0)
  layoutMarginsSubsequentPlot = c(3,0,1,0.08)
  layoutMarginsLastPlot = c(3,0,1,12)
  par(oma = c(1,7,1,7))

  lapply(lpObjects, function(lpObject) {
    solve(lpObject)
    pd <<- poolDifs(gangstaObjects, lpObject, simple = FALSE, byProces = FALSE)
    mt <<- massTransfers(gangstaObjects, lpObject, simple = FALSE, byProcess = FALSE)
    poolOrder <<- abs(seq(-nrow(pd), -1, 1))
    xPosnInit <<- (lpObject %in% lpObjects)
  },
  gangstaObjects
  )


}

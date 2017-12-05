###  Imports gangsta lp file into R
readGangsta.lp = function(lpFile = file.choose()){
  lp.bgc = lpSolveAPI::read.lp(lpFile, type = "lp", verbose = "full")
  return(lp.bgc)
}


### No references to the following two functions (getModelCol.lp and
### showMatrixVals.lp) as of 5 December 2017.  Should we delete this???


### This function allows the user to access the columns of the lp object.
# getModelCol.lp = function(lpObject, colnum) {
#   sapply(1:dim(lpObject)[1], function(i) get.mat(lpObject, i, colnum))
# }

###  This function allows the user to view the values of the lp object.
# showMatrixVals.lp = function(varName, matrix = lpMatrix, lpObject) {
#   varCol = which(colnames(matrix) == varName)
#   varRows = which(getModelCol(lpObject, varCol) != 0)
#   matrix[varRows, varCol]
# }

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

processEnergies =
  function(
    gangstaObjects, lpObject,
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
    rnames = returnDF$procSimple
    returnDF = returnDF[, c("energy", "procType")]
    row.names(returnDF) = rnames
  }
  if(simple) returnDF = subset(returnDF, energy != 0)
  return(returnDF)
}

### Returns the energy from respiration by organism for an lp object
respirationEnergies = function(
  gangstaObjects,
  lpObject){
  nameAttrName = gangstaAttributeName("name")

  respEnergyVarName = gangstaVarName("respEnergy")

  organismClassName = gangstaClassName("org")
  organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
  organismNames = getGangstaAttribute(organisms, nameAttrName)

  respEnergyVarNames = makeGenericVars(organismNames, "respEnergy")

  df.lp = solvedDataFrame.lp(lpObject, simple = F)

  respEnergyVals = df.lp[respEnergyVarNames,]
  names(respEnergyVals) = respEnergyVarNames

  return(respEnergyVals)
}

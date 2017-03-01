## function to pull results of interest from results list
getResultDF = function(results, timestep, DFName) {
  resultDF = results[[timestep]][match(DFName, names(results[[timestep]]))][[1]]
}

## function to aggregate the results of the molTransferVals according to either fromPool or toPool
## returns a vector of the aggregated values
sumTransfersByPool = function(results, timestep, DFName, summedColName, sumFromPools = T) {
  colNames = c("toPool", "fromPool")
  colName = colNames[sumFromPools + 1]

  tableObject = getResultDF(results, timestep, DFName)
  names(tableObject)[names(tableObject) == summedColName] = "col.2.sum"

  summedTransfers = plyr::ddply(tableObject, colName, plyr::summarise, molsTransfered = sum(col.2.sum))

  return(structure(summedTransfers[,"molsTransfered"], names = as.character(summedTransfers[,colName])))
}

getPoolAndTransferVectors = function(results, timestep, poolValDFName, transValDFName, summedColName, poolColName, sumFromPools = T){
  resultsDF = getResultDF(results, timestep, poolValDFName)
  poolNames = row.names(resultsDF)
  poolVals = structure(resultsDF[,poolColName], names = poolNames)
  transferVals = sumTransfersByPool(results, timestep, transValDFName, summedColName, sumFromPools)
  transfersWithZeros = structure(rep(0, length(poolNames)), names = poolNames)
  for(x in names(transferVals)) {
    transfersWithZeros[x] = transferVals[x]
  }
  return(data.frame(poolVals, transfersWithZeros, row.names = poolNames))
}

boxHeights = function(
  results,
  timestep,
  poolValDFName = "poolVals",
  transValDFName = "massTransferVals",
  summedColName = "massTransfered",
  minThickness
) {
  numOfTimestepsInSimulation = length(results)

  resultsList = list()
  if(timestep <= numOfTimestepsInSimulation) {
    resultsList[[1]] =
      getPoolAndTransferVectors(
        results,
        timestep,
        poolValDFName,
        transValDFName,
        summedColName,
        poolColName = "initial",
        sumFromPools = T
      )
  }

  if(timestep >= 2) {
    timestep = timestep - 1
    resultsList[[length(resultsList)+1]] =
      getPoolAndTransferVectors(
        results,
        timestep,
        poolValDFName,
        transValDFName,
        summedColName,
        poolColName = "final",
        sumFromPools = F
      )
  }

  allResultsDF = do.call(cbind, resultsList)

  poolBoxSize = data.frame(size = apply(abs(allResultsDF), 1, max))

  poolBoxSizeNames = row.names(poolBoxSize)[poolBoxSize$size !=0]
  poolBoxSize = data.frame(size = structure(poolBoxSize[poolBoxSize$size !=0,], names = poolBoxSizeNames))

  # ### fudge factor to ensure lines can be seen
  # for(i in 1:length(poolBoxSize$size)){
  #   poolBoxSize$size[i] =
  #     ifelse(poolBoxSize$size[i] < minThickness, minThickness, poolBoxSize$size[i])
  #   }

  # poolNames = gsub("[_]*", "", row.names(poolBoxSize))
  # poolBoxSize$compoundNames = substr(poolNames, 1, nchar(poolNames)-1)
  #
  # boxSize = plyr::ddply(poolBoxSize, "compoundNames", plyr::summarise, boxSize = sum(size))
  #
  # return(boxSize)
  return(poolBoxSize)
}

aggregatedUniqueTransfers = function(
  results,
  timestep,
  DFname = "massTransferVals",
  summedColName = "massTransfered",
  gangstas,
  minThickness = 0.02
) {
  colNames = c("toPool", "fromPool")

  tableObject = getResultDF(results, timestep, DFname)
  names(tableObject)[names(tableObject) == summedColName] = "col.2.sum"

  summedTransfers = plyr::ddply(tableObject, colNames, plyr::summarise, molsTransfered = sum(col.2.sum))
  summedTransfers = summedTransfers[round(summedTransfers$molsTransfered, digits = 12) != 0, ]

  totalMolsTransfered = plyr::ddply(summedTransfers, colNames[2], plyr::summarise, totalMolsTransfered = sum(molsTransfered))
  molsInBoxes = boxHeights(results, timestep)
  zeroMols = row.names(molsInBoxes)[!(row.names(molsInBoxes) %in% totalMolsTransfered$fromPool)]
  if(length(zeroMols) > 0) {
    zeroMolsTransfered =
      data.frame(
        fromPool =  zeroMols,
        totalMolsTransfered = 0
      )
    totalMolsTransfered = rbind(totalMolsTransfered, zeroMolsTransfered)
  }
  molsInBoxes =
    molsInBoxes[with(molsInBoxes, match(totalMolsTransfered$fromPool, row.names(molsInBoxes))),]
  names(molsInBoxes) = totalMolsTransfered$fromPool
  molsRemaining = molsInBoxes - totalMolsTransfered$totalMolsTransfered

  toPool = rep(NA, length(molsRemaining))
  fromPool = rep(NA, length(molsRemaining))
  molsTransfered = rep(-999, length(molsRemaining))
  for(i in 1:length(molsRemaining)){
    toPool[i] = names(molsRemaining)[i]
    fromPool[i] = names(molsRemaining)[i]
    molsTransfered[i] = molsRemaining[i]
  }
  newRowsToAdd = data.frame(toPool = toPool, fromPool = fromPool, molsTransfered = molsTransfered)
  newRowsToAdd = newRowsToAdd[round(newRowsToAdd$molsTransfered, digits = 12) != 0, ]

  summedTransfers = rbind(summedTransfers, newRowsToAdd)

  compounds = subsetGangstas(gangstas, "class", gangstaClassName("comp"))
  InfinteCompoundLogicalVector = getGangstaAttribute(compounds, gangstaAttributeName("InfinteCompound"))
  compoundNames = getGangstaAttribute(compounds, gangstaAttributeName("name"))
  names(InfinteCompoundLogicalVector) = compoundNames
  InfinteCompounds = compoundNames[InfinteCompoundLogicalVector]

  fromNames = as.character(summedTransfers$fromPool)
  toNames = as.character(summedTransfers$toPool)
  fromToIdentical = (fromNames == toNames)
  compoundNamesInTransfersDF = substr(fromNames, 1, nchar(fromNames)-2)
  isSourceSink = (compoundNamesInTransfersDF %in% InfinteCompounds)

  InfinteCompoundRemovalCriteria = (fromToIdentical & isSourceSink)
  keep = !(InfinteCompoundRemovalCriteria)
  summedTransfers = summedTransfers[keep, ]

  return(summedTransfers)
}

simplifyDataForAlluvialPlot = function(results, gangstas) {
  numOfTimesteps = length(results)
  # boxSizes = lapply(1:(numOfTimesteps+1), function(ts) boxHeights(results, ts))
  molTransfers = lapply(1:numOfTimesteps, function(ts) aggregatedUniqueTransfers(results, ts, gangstas = gangstas))
  # return(list(boxes = boxSizes, arrows = molTransfers))
  return(list(arrows = molTransfers))

}

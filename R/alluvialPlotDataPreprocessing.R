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

  # Get results from results list for the time step indicated
  tableObject = getResultDF(results, timestep, DFname)
  # Rename the "massTransfered" column name
  names(tableObject)[names(tableObject) == summedColName] = "col.2.sum"

  # Combine the transfers where the toPools and fromPools are the same
  summedTransfers = plyr::ddply(tableObject, colNames, plyr::summarise, molsTransfered = sum(col.2.sum))

  # Sum the mols transfered by fromPool
  totalMolsTransfered = plyr::ddply(summedTransfers, colNames[2], plyr::summarise, totalMolsTransfered = sum(molsTransfered))
  # Get the heights of the boxes for each pool
  molsInBoxes = boxHeights(results, timestep)

  # Get the fromPools not in the molsInBoxes data frame and add them to the molsInBoxes data frame
  poolsWithZeroMolsInBoxes = totalMolsTransfered$fromPool[!(totalMolsTransfered$fromPool %in% row.names(molsInBoxes))]
  if(length(poolsWithZeroMolsInBoxes) > 0 ) {
    zeroMolsInBoxes = as.data.frame(rep(0, length(poolsWithZeroMolsInBoxes)))
    row.names(zeroMolsInBoxes) = poolsWithZeroMolsInBoxes
    names(zeroMolsInBoxes) = "size"
    molsInBoxes = rbind(molsInBoxes, zeroMolsInBoxes)
  }
  molsInBoxes$size = round(molsInBoxes$size, digits = 10)

  # Get the molsInBoxes that are not listed in the fromPools and add them to the totalMolsTransfered data frame
  poolsWithZeroTransfers = row.names(molsInBoxes)[!(row.names(molsInBoxes) %in% totalMolsTransfered$fromPool)]
  if(length(poolsWithZeroTransfers) > 0) {
    zeroTransfers = data.frame(
      fromPool = poolsWithZeroTransfers,
      totalMolsTransfered = rep(0, length(poolsWithZeroTransfers))
      )
    totalMolsTransfered = rbind(totalMolsTransfered, zeroTransfers)
  }
  totalMolsTransfered$totalMolsTransfered = round(totalMolsTransfered$totalMolsTransfered, digits = 10)

  # order the mols in boxes data frame by pool name
  molsInBoxesVect = molsInBoxes$size[order(row.names(molsInBoxes))]
  names(molsInBoxesVect) = row.names(molsInBoxes)[order(row.names(molsInBoxes))]

  # order the totalMolsTransfered by fromPool name
  totalMolsTrans = totalMolsTransfered$totalMolsTransfered[order(as.character(totalMolsTransfered$fromPool))]
  names(totalMolsTrans) = totalMolsTransfered$fromPool[order(as.character(totalMolsTransfered$fromPool))]

  # print(timestep)

  if(!identical(names(molsInBoxesVect),names(totalMolsTrans))) {
    stop("The molsInBoxesVect doesn't match the totalMolsTrans!")
    }

  # Calculate the number of mols remaining in boxes at the end of the time step
  molsRemaining = round(molsInBoxesVect - totalMolsTrans, digits = 10)

  # Add molsRemaining to the summedTransfers data fram where the fromPool and
  # toPool are identical.  Essentially, this allows flow through between
  # identical pools.
  newRowsToAdd =
    data.frame(
      toPool = names(molsRemaining),
      fromPool = names(molsRemaining),
      molsTransfered = molsRemaining
      )
  summedTransfers = rbind(summedTransfers, newRowsToAdd)

  # This next block of code gets rid of flow through between time steps for
  # Infinite Compounds.

  # Identify the infinite compounds
  compounds = subsetGangstas(gangstas, "class", gangstaClassName("comp"))
  InfiniteCompoundLogicalVector = getGangstaAttribute(compounds, gangstaAttributeName("InfiniteCompound"))
  compoundNames = getGangstaAttribute(compounds, gangstaAttributeName("name"))
  names(InfiniteCompoundLogicalVector) = compoundNames
  InfiniteCompounds = compoundNames[InfiniteCompoundLogicalVector]

  fromNames = as.character(summedTransfers$fromPool)
  toNames = as.character(summedTransfers$toPool)
  fromToIdentical = (fromNames == toNames)
  compoundNamesInTransfersDF = substr(fromNames, 1, nchar(fromNames)-2)
  isSourceSink = (compoundNamesInTransfersDF %in% InfiniteCompounds)

  # Get rid of flow through for infinite compounds
  InfiniteCompoundRemovalCriteria = (fromToIdentical & isSourceSink)
  keep = !(InfiniteCompoundRemovalCriteria)
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

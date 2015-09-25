###  Imports gangsta lp file into R
readGangsta.lp = function(lpFile = file.choose()){
  lp.bgc = lpSolveAPI::read.lp(lpFile, type = "lp", verbose = "full")
  return(lp.bgc)
}

### This function allows the user to access the columns of the lp object.
getModelCol.lp = function(gangsta.lp, colnum) {
  sapply(1:dim(lp.model)[1], function(i) get.mat(lp.model, i, colnum))
}

###  This function allows the user to view the values of the lp object.
showMatrixVals.lp = function(varName, matrix = lpMatrix, lpObject = gangsta.lp) {
  varCol = which(colnames(matrix) == varName)
  varRows = which(getModelCol(lpObject, varCol) != 0)
  matrix[varRows, varCol]
}

###  This block of code builds a matrix of decision variables and constraints
viewMatrix.lp = function(gangsta.lp){
  constraintVarNames = dimnames(gangsta.lp)[[1]]
  decisionVarNames = dimnames(gangsta.lp)[[2]]
  nvars = length(decisionVarNames)
  nconstraints = length(constraintVarNames)
  lpMatrix = matrix(mapply(function(i,j) lpSolveAPI::get.mat(gangsta.lp, i, j),
                           i=1:nconstraints, j=rep(1:nvars, each=nconstraints)),
                    nrow = nconstraints)
  colnames(lpMatrix) = decisionVarNames
  rownames(lpMatrix) = constraintVarNames
  return(lpMatrix)
}

### View output in data frame
solvedDataFrame.lp = function(gangsta.lp, simple = TRUE) {
  variableOutput = lpSolveAPI::get.variables(gangsta.lp)[order(dimnames(gangsta.lp)[[2]])]
  names(variableOutput) = dimnames(gangsta.lp)[[2]][order(dimnames(gangsta.lp)[[2]])]
  variableOutput = data.frame(variableOutput)
  if(simple == FALSE){
    return(variableOutput)
  } else {
    variableOutputSimple = subset(variableOutput, variableOutput != 0)
    variableOutputSimple[order(row.names(variableOutputSimple)),]
    return(variableOutputSimple)
  }
}


poolDifs = function(gangstaObjects, lp, simple = T) {
  pools = subsetGangstas(gangstaObjects, "class", gangstaClassName("pool"))
  poolNames = getGangstaAttribute(pools, gangstaAttributeName("name"))

  initPoolNames = paste0(poolNames, ".", gangstaVarName("startSuffix"))
  finalPoolNames = paste0(poolNames, ".", gangstaVarName("endSuffix"))

  df.lp = solvedDataFrame.lp(lp, simple = F)

  startVals = df.lp[initPoolNames,]
  endVals = df.lp[finalPoolNames,]
  allDF = data.frame(change = endVals - startVals, initial = startVals, final = endVals)
  row.names(allDF) = poolNames
  if(simple) allDF = subset(allDF, (allDF$initial != 0) | (allDF$final != 0))
#  col.names(appDF) = c("initial", "final", "change")
  return(allDF)
}

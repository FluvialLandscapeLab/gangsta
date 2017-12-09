#User provides:
#
dummyf = function() {
  return(1)
}

slopesToUpdate =
 list(
   list(constraint = c("Het.finalMolecules", "Het.respirationEnergy"), funct = dummyf, initValue = NULL),
   list(constraint = c("Aut.finalMolecules", "Aut.respirationEnergy"), funct = dummyf, initValue = NULL),
   list(constraint = c("Met.finalMolecules", "Met.respirationEnergy"), funct = dummyf, initValue = NULL)
 )

# constraint is identified as the row in the lpmatrix where all "vars" have a
# non-zero slope, but other vars have a zero slope slope associated with "var1"
# is always replaced by value from function.

# variablesToUpdate =
#   list(
#     list(variable = "variableName", function = f, initValue = NULL),
#     list(variable = "variableName", function = f, initValue = NULL)
#   )

# if initValue is numeric(0), rely on initialized value within the lp code
# if initValue is NULL, run the function to initialize the value
# if initValue is a numeric value, set the variable to the value.

# drivingValues = dataframe(drivingVar1 = c(...), drivingVar2 = c(...), rownames = 0:nsteps)

GANGSTARap = function(lpModel, slopesToUpdate, variablesToUpdate, drivingValues) {

  # helper function that extracts named values from sublist of lists.
  getListParts = function(lst, nm) lapply(lst, "[[", nm)

  getConstraintRow = function(constraintID) {
    badName = !constraintID %in% colnames(lpModelMatrix)
    if(any(badName)) {
      stop("Variable(s) named '", paste(constraintID[badName], sep = "', '"), "' are not in the model.")
    }
    keyColumns = match(constraintID, colnames(lpModelMatrix))
    hasRequiredSlopes = apply(lpModelMatrix, 1, function(x) !any(x[keyColumns] == 0))
    hasRequiredZeros = apply(lpModelMatrix, 1, function(x) !any(x[-keyColumns] != 0))
    constraintRow = which(hasRequiredZeros & hasRequiredSlopes)
    if(length(constraintRow) == 0) stop("No constraint found that uses only the variable(s): '", paste(constraintID, sep="', '"), "'")
    if(length(constraintRow) > 1) stop("More than one constraint found that uses only the variable(s): '", paste(constraintID, sep="', '"), "'")
    return(constraintRow)
  }

  constraintIDs = getListParts(slopesToUpdate, "constraint")
  slopeFunctions = getListParts(slopesToUpdate, "funct")
  initvalues = getListParts(slopesToUpdate, "initValues")

  lpModelMatrix = viewMatrix.lp(lpModel)

  constraintRows = sapply(constraintIDs, getConstraintRow)
  constrainColumns = match(sapply(constraintIDs, "[", 1), colnames(lpModelMatrix))
  newSlopes = data.frame(row = constraintRows, column = constrainColumns)
  newSlopes$funct = slopeFunctions
  return(newSlopes)

#  this line will now update all slopes with whatever function is specified to calculate the slope.
#  mapply(function(r, c, f) set.mat(mylp, r, c, do.call(f, args = list())), test$row, test$column, test$funct)
}

#newSlopes data.frame with lpmatrix row, column, and function
updateSlope = function(lpModel, newSlopes, rapperEnvir = parent.frame()) {
  mapply(
    set.mat,
    i = newSlopes$row,
    j = newSlopes$column,
    value = do.call(newSlopes$funct, args = list(), envir = rapperEnvir),
    moreArgs = list(lprec = lpModel)
  )
}

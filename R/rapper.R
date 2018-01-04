#User provides:
# slopesToUpdate =
#  list(
#    list(constraint = c("Var1", "Var2", "VarN"), funct = f, initValue = NULL)
#  )
# constraint is identified as the row in the lpmatrix where all "vars" have a
# non-zero slope, but other vars have a zero slope. The slope associated with "var1"
# is always replaced by value from function.


# parametersToUpdate =
#   list(
#     list(parameter = "parameterName", funct = f, initValue = NULL),
#     list(parameter = "parameterName", funct = f, initValue = NULL)
#   )
# parametersToUpdate consists of any value that is set equal to a constant with upper
# and lower constraints equal to that value. (i.e. Het.initialMolecules = 0)

# if initValue is numeric(0), rely on initialized value within the lp code
# if initValue is NULL, run the function to initialize the value
# if initValue is a numeric value, set the parameter to the value.

# drivingValues = dataframe(drivingVar1 = c(...), drivingVar2 = c(...), rownames = 0:nsteps)

# EXAMPLE
# lpModel = readGangsta.lp()
#  dummyf = function() {
#    return(1)
#  }
#
#  slopesToUpdate =
#    list(
#     list(constraint = c("Het.finalMolecules", "Het.respirationEnergy"), funct = dummyf, initValue = NULL),
#     list(constraint = c("Aut.finalMolecules", "Aut.respirationEnergy"), funct = dummyf, initValue = NULL),
#     list(constraint = c("Met.finalMolecules", "Met.respirationEnergy"), funct = dummyf, initValue = NULL)
#   )
#
#
#  parametersToUpdate =
#    list(
#      list(parameter = "Het.initialMolecules", funct = dummyf, initValue = NULL),
#      list(parameter = "Aut.initialMolecules", funct = dummyf, initValue = NULL)
#    )


# GangstaRap is a function that builds a rapper around GANGSTA which allows the GANGSTA environment
# to recallibrate to changes following each timestep. The values that need to change are defined
# in slopesToUpdate and parametersToUpdate. The values are changed according to the function defined in
# slopesToUpdate and parametersToUpdate (funct = f). The drivingValues is a dataframe of values to be input
# into function "f". The lpModel is the lp file created by the initiation of the GANGTSA model. Changes
# made by the GANGSTARap function will be made directly to the lp file specified.
GANGSTARap = function(lpModel, slopesToUpdate, parametersToUpdate, drivingValues, envir = parent.frame()) {

  # helper function that extracts named values from sublist of lists.
  getListParts = function(lst, nm) lapply(lst, "[[", nm)
  # helper function that identifies if any variables passed through slopeToUpdate or constraintsToUpdate
  # are not in the model and ends the GANGSTARap function if this is true for a variable name.
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
  updateSlopes(lpModel, newSlopes, envir = envir)

# this line will now update all slopes with whatever function is specified to calculate the slope.
# mapply(function(r, c, f) set.mat(mylp, r, c, do.call(f, args = list())), test$row, test$column, test$funct)

# gets list parts from parametersToUpdate
  parameterIDs =getListParts(parametersToUpdate, "parameter")
  parameterFunctions = getListParts(parametersToUpdate, "funct")
  parameterInitVals = getListParts(parametersToUpdate, "initValue")

# Gets location of parameters from lpModel for parametersToUpdate
  parameterIndexes = (match(parameterIDs, dimnames(lpModel)[[2]]))
# Gets existing bounds for parametersToUpdate. If bounds do not match, returns an error.
  existingParameterBounds = lpSolveAPI::get.bounds(lpModel, columns = parameterIndexes)
  if(!(identical(existingParameterBounds$lower, existingParameterBounds$upper))) {
    stop("Upper and lower bounds for parameter bounds in existing lp models don't match.")
  }
  newParameters = data.frame(row = parameterIndexes)
  newParameters$funct = parameterFunctions
  updateParameters(lpModel, newParameters, envir = envir)


}
# end of GangstaRap

# set new lp.matrix values using results from a function
set.matWithFunction = function(lpModel, i, j, fun, args = list(), envir = parent.frame()){
  set.mat(
    lpModel,
    i = i,
    j = j,
    value = do.call(fun, args = args, envir = envir)
  )
}

#update newSlopes data.frame with lpmatrix row, column, and function
updateSlopes = function(lpModel, newSlopes, envir = parent.frame()) {
  mapply(
    set.matWithFunction,
    i = newSlopes$row,
    j = newSlopes$column,
    fun = newSlopes$funct,
    MoreArgs = list(lpModel= lpModel, envir = envir)
  )
}

#set new lp.bounds using results from a function
set.boundsWithFunction = function(lpModel, column, fun, upperFun = fun, envir = parent.frame()){
  set.bounds(
    lpModel,
    lower = do.call(fun, args = list(), envir = envir),
    upper = do.call(upperFun, args = list(), envir = envir),
    columns = column
  )
}

# update parameters with
updateParameters = function(lpModel, newParameters, envir = parent.frame()){
  mapply(
    set.boundsWithFunction,
    column = newParameters$row,
    fun = newParameters$funct,
    MoreArgs = list(lpModel = lpModel, envir = envir)
  )
}





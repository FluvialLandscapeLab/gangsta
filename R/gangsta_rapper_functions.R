getValueNamesToUpdate = function(lpModel){
  # get variable names from model
  modelVarNames = dimnames(lpModel)[[2]]
  
  # make a list of all of the names that need to be updated between model runs by
  # finding all names that have ".intialMolecules" in them.
  initialValueNames = modelVarNames[grep(".initialMolecules", modelVarNames)]
  # Ooops!  Except remove the Ox. and Hx. variables
  initialValueNames = initialValueNames[-grep("x.", initialValueNames)]
  
  # create the names of that variables that have the values to be used in updating
  # "updateName" variables.  (i.e. the ".finalMolecules" variables have the values
  # from the prior time step.)
  finalValueNames = sub("initial", "final", initialValueNames)
  return(list(initialValueNames = initialValueNames, finalValueNames = finalValueNames))
}

getLeakInValueNames = function(initialValueNames, drivingValues) {
  addToValueNames = paste0("add.to.", initialValueNames)
  hasLeakIn = addToValueNames %in% names(drivingValues)
  addToValueNames =
    structure(
      addToValueNames[hasLeakIn],
      names = initialValueNames[hasLeakIn]
    )
  return(addToValueNames)
}

get.matRow <- function(rowNum, lpModel){
  sapply(1:ncol(lpModel), function(colNum) {lpSolveAPI::get.mat(i = rowNum, j = colNum, lprec = lpModel)})
}

# helper function that returns the row number of a constraint.  The
# "constraintID" is a vector of variables names used by the constraint.
findConstraintRowAndColumn = function(constraintID, lpModel) {
  # check to be sure all variable names in the constraintID are in the model.
  badName = !constraintID %in% dimnames(lpModel)[[2]]
  if(any(badName)) {
    stop("Variable(s) named '", paste(constraintID[badName], collapse = "', '"), "' are not in the model.")
  }
  ## Libby changed this line of code, 9-23-18
  #  modelVarNames = get.lpModelVarNames(lpModel)
  modelVarNames = dimnames(lpModel)[[2]]
  # the model has a matrix where each row represents a constraint and each
  # column represents the slopes associated with the variable for the
  # constraint. So, identifying a constraint row means that we find the row
  # with no-zero slopes for each variable in the constraintID and zero slope
  # for every variable not in the constraintID.
  
  # the follow columns need to have non-zero slopes
  constraintColumns = match(constraintID, modelVarNames)
  constraintColumn = constraintColumns[1]
  
  # return row number matrix row where constraintColumns are non-zero and other columns are zero in a constraint row
  # What is get.matRow? Added a new function below. (Libby, 9-23-18)
  constraintRow =
    which(
      sapply(
        1:nrow(lpModel),
        function(rowNum) {
          nonZeroColumns <- which(get.matRow(rowNum, lpModel) != 0.0)
          return(identical(sort(as.integer(nonZeroColumns)), sort(as.integer(constraintColumns))))
        }
      )
    )
  # makes sure there is one and only one constraint row that matches constraintID
  if(length(constraintRow) == 0) stop("No constraint found that uses only the variable(s): '", paste(constraintID, sep="', '"), "'")
  if(length(constraintRow) > 1) stop("More than one constraint found that uses only the variable(s): '", paste(constraintID, sep="', '"), "'")
  return(structure(c(constraintRow, constraintColumn), names = c("i", "j")))
}

# getModelMatrix returns a matrix for the solved lpModel
getModelMatrix = function(envir){
  mat = matrix(
    mapply(
      lpSolveAPI::get.mat,
      i = rep(1:nrow(envir$lpModel), ncol(envir$lpModel)),
      j = rep(1:ncol(envir$lpModel), each = nrow(envir$lpModel)),
      MoreArgs = list(lprec = envir$lpModel)
    ),
    nrow = nrow(envir$lpModel),
    dimnames = dimnames(envir$lpModel)
  )
}

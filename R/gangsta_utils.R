#' gangsta utils
#'
#' \code{writeGangstaModel} instantiates a \code{gangsta}-derived model into
#' computer code.
#'
#' @param gangstaObjects \code{gangstaObjects} are the list of objects of class
#'   \code{gangsta} created with functions such as \code{compoundFactory} and
#'   \code{processFactory}.
#' @param file The default argument for \code{file} is \code{file.choose()};
#'   using the default will enable the end user to overwrite an existing file or
#'   create a new file if one does not already exist.  However, the end user can
#'   also specify a file on their local drive. In either case, files generated
#'   using \code{writeGangstaModel} should be saved with an ".lp" extension.
#'
#' @return \code{writeGangstaModel} returns a file with the simulation model
#'   code, formatted for use in lpSolve (http://lpsolve.sourceforge.net/5.5/).
#' @export
writeGangstaModel = function(gangstaObjects, file = file.choose()) {
  expressions = makeExpressions(gangstaObjects) # create the expressions
  expressions = formatExpressions(expressions) # add the semicolon to the expressions
  file.create(file) # create the .lp file
  write(expressions, file) # write the .lp file
}

makePoolNames = function(compoundNames, elementNames = names(compoundNames), gangstaObjects = NULL) {
  poolNames = paste0(compoundNames, "_", elementNames)
  if(!is.null(gangstaObjects)) {
    gangstasExist(gangstaObjects, poolNames, gangstaClassName("pool"))
  }
  return(poolNames)
}

# Return the names of all of the paramers of the processFactory() function that
# don't have default values.  Test for new default value is whether or not
# eval() throws an error.  If so, there is no default.
processSpecRequiredNames = function() {
  argHasNoDefault = function(x) {
    evaluated = tryCatch(
      eval(x),
      error = function(e) return(">>NODEF<<")
    )
    if(identical(evaluated, ">>NODEF<<")) {
      return(T)
    } else {
      return(F)
    }
  }

  processFactoryParams = formals(processFactory)
  processFactoryParams = processFactoryParams[2:length(processFactoryParams)]
  hasNoDefault = sapply(processFactoryParams, argHasNoDefault)
  names(hasNoDefault) = names(processFactoryParams)
  return(hasNoDefault)
}

checkProcessSpecNames = function(processSpecNames, additionalValidNames = NULL) {
  namesAreRequired = processSpecRequiredNames()
  validNames = c(names(namesAreRequired), additionalValidNames)
  isInvalid = !(processSpecNames %in% validNames)
  if(any(isInvalid)) {
    stop("The following names are not permitted in a process specification: ", paste0(processSpecNames[isInvalid], collapse = "; "), "\n  (Valid names are: ", paste0(validNames, collapse = "; "), ")")
  }

  # hint: namesAreRequired is a boolean vector, so this next line returnes the
  # names where the vector value is T
  requiredNames = names(namesAreRequired[namesAreRequired])
  isMissing = !(requiredNames %in% processSpecNames)
  if(any(isMissing)) {
    stop("The following names are required in a process specification but are missing: ", paste0(requiredNames[isMissing], collapse = "; "))
  }

  return(T)
}




fromToPair = function(gangstaObjects, fromName, toName) {
  targetPoolNames = c(fromName, toName)
  gangstasExist(gangstaObjects, targetPoolNames, "pool")
  targetPools = getGangstas(gangstaObjects, targetPoolNames)
  if(identical(targetPools[[1]], targetPools[[2]])) {
    stop(paste0("A requested transfer has '", fromName, "' as both the 'to' and 'from' pools.  Tranformations can't connect a pool to itself."))
  }
  elements = getGangstaAttribute(targetPools, "elementName")
  if(length(unique(elements))>1){
    stop("A transfer can only link pools that contain the same element.")
  }
  names(targetPools) = c("from", "to")
  return(targetPools)
}

gangstaVarName = function(varTag) {
  return(getOption("gangsta.vars")[varTag])
}

gangstaClassName = function(classTag) {
  return(getOption("gangsta.classes")[classTag])
}

gangstaAttributeName = function(attributeTag) {
  return(getOption("gangsta.attributes")[attributeTag])
}

gangstaVarTags = function() {
  return(names(getOption("gangsta.vars")))
}

gangstaClassTags = function() {
  return(names(getOption("gangsta.classes")))
}

gangstaAttributeTags = function() {
  return(names(getOption("gangsta.attributes")))
}



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
  #  targetPools = lapply(targetPoolNames, subsetGangstas, gangstaObjects = pools, attribName = "name")
  targetPools = getGangstas(gangstaObjects, targetPoolNames)
  if(identical(targetPools[[1]], targetPools[[2]])) {
    stop(paste0("A requested transformation has '", fromName, "' as both the 'to' and 'from' pools.  Tranformations can't connect a pool to itself."))
  }
  elements = getGangstaAttribute(targetPools, "elementName")
  if(length(unique(elements))>1){
    stop("A transformation can only link pools that contain the same element.")
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

writeGangstaModel = function(equations, file = file.choose()) {
  file.create(file)
  equations = formatEquations(equations)
  write(equations, file)
}

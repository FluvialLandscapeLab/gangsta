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

expandMultiprocessSpec = function(multiprocessSpec) {
  checkProcessSpecNames(names(multiprocessSpec), "processSuffix")

  # if a "processSuffix" is passed in multiprocessSpec, we we don't want to
  # return it from this function. Realizing the T = 1, and F = 0, calculated the
  # number of elements to return.
  nListElementsToReturn = length(multiprocessSpec) - !is.null(multiprocessSpec[["processSuffix"]])

  # add "" for organismName and processSuffix if not already in list
  defaults = list(organismName = "", processSuffix = "")
  notInList = !(names(defaults) %in% names(multiprocessSpec))
  spec = c(multiprocessSpec, defaults[notInList])

  # move processSuffix to end of list
  suffixList = spec["processSuffix"]
  spec["processSuffix"] = NULL
  spec = c(spec, suffixList)

  # use names of spec to create variables with associated values
  for(varName in names(spec)) {
    assign(varName, spec[[varName]])
  }

  # expand the elements of spec
  organismName = rep(organismName, each = length(processSuffix))
  name = paste0(organismName, name, processSuffix)
  if(!is.null(spec[["limitToInitMols"]])) limitToInitMols = rep(limitToInitMols, each = length(processSuffix))
  fromCompoundNames = do.call(mapply, c(list(FUN = "c", SIMPLIFY = F), fromCompoundNames))
  toCompoundNames = do.call(mapply, c(list(FUN = "c", SIMPLIFY = F), toCompoundNames))
  molarTerms = do.call(mapply, c(list(FUN = "c", SIMPLIFY = F), molarTerms))

  # Store the values of the expanded elements in a list
  paramList = lapply(names(spec)[1:nListElementsToReturn], function(x) eval(parse(text = x)))
  # Now use do.call to call mapply on the "list" fuction.  This is the same as
  # executing: mapply(FUN = list, name, energyTerm, fromCompoundNames,
  # toCompoundNames, molarTerms, organismName, limitToInitMols)
  # which obviously(?) creates multiple processSpecs from the recyled
  # combinations of the expanded elements of the multiProcessSpec
  expandedSpecs = do.call(mapply, c(list(FUN = "list", SIMPLIFY = F), paramList))
  expandedSpecs = lapply(
    expandedSpecs,
    function(x) {
      names(x) = names(spec)[1:nListElementsToReturn]
      return(x)
    }
  )
  expandedSpecs = lapply(expandedSpecs, processSpec)
  names(expandedSpecs) = name
  return(expandedSpecs)
}


replaceNAWithMolarRatio = function(molarTerms, fromPoolNames, fromCompoundNames, gangstaObjects) {

  isNA = is.na(molarTerms)
  if(any(isNA)) {
    fromCompounds = getGangstas(gangstaObjects, fromCompoundNames)
    refPoolNames = getGangstaAttribute(fromCompounds, gangstaAttributeName("refPool"))
    refPoolIndexes = match(refPoolNames, fromPoolNames)
    missingRefPools = (is.na(refPoolIndexes) & isNA)
    if(any(missingRefPools)) {
      stop("When molarTerms are 'NA' in a process specification, the molarTerm of the reference pool must be specified in the vector.\n  The molarTerm(s) for ",
           paste0(fromPoolNames[match(refPoolNames[missingRefPools], refPoolNames)], collapse = " and "),
           " = NA, but expected molarTerm(s) for reference pool(s) ",
           paste0(refPoolNames[missingRefPools], collapse = " and "),
           " are not in the molarTerms vector.\n  Error occurred where fromPools = ",
           paste0(fromPoolNames, collapse = "; ")
      )
    }
    # in case the refPool is specified as the fromPool for more than one
    # transformation (e.g., DOM -> CO2 and DOM -> CH4 in denitrification), we
    # sum the molarTerms for all instances of each referencePool in fromPools
    refPoolMolarTerms = sapply(refPoolNames, function(rP) (sum(molarTerms[rP == fromPoolNames])))
    isNARefMolarTerms = (is.na(refPoolMolarTerms) & isNA)
    if(any(isNARefMolarTerms)) {
      stop("The molarTerms for reference pools can not be 'NA' in a process specification.\n  The molarTerm(s) for ",
           paste0(unique(refPoolNames[isNARefMolarTerms]), collapse = " and "),
           " were specified as 'NA'.\n  Error occurred where fromPools = ",
           paste0(fromPoolNames, collapse = "; ")
      )
    }
    fromPools = getGangstas(gangstaObjects, fromPoolNames)
    fromPoolMolarRatioList = getGangstaAttribute(fromPools, gangstaAttributeName("molRatio"))
    fromPoolMolarRatios = sapply(
      fromPoolMolarRatioList,
      function(x) {
        if (is.null(x)) return(NA)
        return(x)
      }
    )
    molarTerms[isNA] = refPoolMolarTerms[isNA] * fromPoolMolarRatios[isNA]
  }
  return(molarTerms)
}

replaceDotWithOrganism = function(compoundNames, organismName) {
  compoundNames[compoundNames == "."] = organismName
  return(compoundNames)
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
  equations = sapply(equations, paste0, ";")
  write(equations, file)
}

processSpec = function(listOfProcessFactoryArgs) {
  checkProcessSpecNames(names(listOfProcessFactoryArgs))
  class(listOfProcessFactoryArgs) = c("processSpec", "gangsta")
  return(listOfProcessFactoryArgs)
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

replaceDotWithOrganism = function(compoundNames, organismName) {
  compoundNames[compoundNames == "."] = organismName
  return(compoundNames)
}





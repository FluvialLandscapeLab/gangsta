processSpec = function(listOfProcessFactoryArgs) {
  checkProcessSpecNames(names(listOfProcessFactoryArgs))
  class(listOfProcessFactoryArgs) = c("processSpec", "gangsta")
  return(listOfProcessFactoryArgs)
}

enList = function(x) {
  if(!is.list(x)) x = list(x)
  return(x)
}
expandMultiprocessSpec = function(multiprocessSpec) {
  checkProcessSpecNames(names(multiprocessSpec), "processSuffix")

  # if a "processSuffix" is passed in multiprocessSpec, we we don't want to
  # return it from this function. Realizing the T = 1, and F = 0, calculate the
  # number of list elements to return.
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

  # ------expand the elements of spec
  # 1) replicate organism names for each proccessSuffix
  organismName = rep(organismName, each = length(processSuffix))
  name = paste0(organismName, name, processSuffix)
  # 2) if limitToMols is specified, replicate that as well
  if(!is.null(spec[["limitToInitMols"]])) limitToInitMols = rep(limitToInitMols, each = length(processSuffix))
  # 3) create a list of from vectors where the length of the list is equal to
  # the highest number of specified fromCompounds and the other fromCompounds
  # are recycled
  fromCompoundNames = lapply(fromCompoundNames, enList)
  fromCompoundNames = do.call(mapply, c(list(FUN = "list", SIMPLIFY = F), fromCompoundNames))
  # 4) same for toCompounds and molar terms
  toCompoundNames = lapply(toCompoundNames, enList)
  toCompoundNames = do.call(mapply, c(list(FUN = "list", SIMPLIFY = F), toCompoundNames))
  molarTerms = do.call(mapply, c(list(FUN = "c", SIMPLIFY = F), molarTerms))
  # ------end exansion

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
  expandedSpecs = lapply(expandedSpecs, expandTransferGroups)
  expandedSpecs = lapply(expandedSpecs, processSpec)
  names(expandedSpecs) = name

  for (i in 1:length(expandedSpecs)) {
    expandedSpecs[[i]]$fromCompoundNames = replaceDotWithOrganism(expandedSpecs[[i]]$fromCompoundNames, expandedSpecs[[i]]$organismName)
    expandedSpecs[[i]]$toCompoundNames = replaceDotWithOrganism(expandedSpecs[[i]]$toCompoundNames, expandedSpecs[[i]]$organismName)
  }

  return(expandedSpecs)
}

expandTransferGroups = function(specList) {
  # return the pairs of from and to compounds by matching from and to vectors
  # and recycling when vector lengths are different
  transferOptions =
    mapply(
      function(f, t) mapply(c, f, t, SIMPLIFY = F),
      specList$fromCompoundNames,
      specList$toCompoundNames,
      SIMPLIFY = F
    )
  # restructure pairs in the list of lists to an array; first column is the
  # expandedFroms, second is expandedTos
  transferOptionMatrix = do.call(rbind, unlist(transferOptions, recursive = F))
  expandedFrom = transferOptionMatrix[,1]
  expandedTo = transferOptionMatrix[,2]

   # calculate number of times to rep each from, to, and molarTerm
  repTimes = pmax(sapply(specList$fromCompoundNames, length), sapply(specList$toCompoundNames, length))
  # rep the names and assign to the expanded Froms and Tos
  names(expandedFrom) = rep.int(names(specList$fromCompoundNames), repTimes)
  names(expandedTo) = rep.int(names(specList$toCompoundNames), repTimes)

  # expand the molar terms based on the length of each of the from and to
  # compound vectors
  expandedMolarTerms = rep.int(specList$molarTerms, times = repTimes)
  names(expandedMolarTerms) = rep.int(names(specList$molarTerms), times = repTimes)

  # calculate the transferOption indexes.  Don't try to figure this out.
  nTrans = length(expandedFrom)
  startIdx = cumsum(c(1, repTimes[-length(repTimes)]))
  endIdx = cumsum(repTimes)
  idxList = mapply(":", startIdx, endIdx, USE.NAMES = F, SIMPLIFY = F)
  names(idxList) = names(specList$molarTerms)

  # overwrite original froms and tos in the specList
  specList$fromCompoundNames = expandedFrom
  specList$toCompoundNames = expandedTo
  specList$molarTerms = expandedMolarTerms

  #add the transferOptions to the specList
  specList$transferOptions = idxList

  # # expand the molarTerms to be 1:1 with number of transfers
  # expandedMolarTerms = rep.int(specList$molarTerms, times = sapply(idxList, length))
  # names(expandedMolarTerms) = rep.int(names(specList$molarTerms), times = sapply(idxList, length))
  #
  # # make each molarTerm a list containing the value and the transferGroup indexes
  # newMolarTerms = list(groupValues = specList$molarTerms, groupIdx = idxList, expandedValues = expandedMolarTerms)
  #
  # # newMolarTerms = mapply(list, specList$molarTerms, idxList, SIMPLIFY = F)
  # # name the molarTerm elements "value" and "groupIdx"
  # # newMolarTerms =
  # #   lapply(
  # #     newMolarTerms,
  # #     function(l, n) {
  # #       names(l) = n
  # #       return(l)
  # #     },
  # #     n = c("value", "groupIdx")
  # #   )
  #
  # # update specList
  # specList$molarTerms = newMolarTerms
  return(specList)
}

replaceDotWithOrganism = function(compoundNames, organismName) {
  compoundNames[compoundNames == "."] = organismName
  return(compoundNames)
}





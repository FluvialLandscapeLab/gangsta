makePoolName = function(compoundName, elementName) {
  paste0(compoundName, "_", elementName)
}

replaceDotWithOrganism = function(compoundNames, organismName) {
  compoundNames = lapply(
    compoundNames,
    function(c) {
      c[c == "."] = organismName
      return(c)
    }
  )
  return(compoundNames)
}

## requires a named list; elementNames are the names of the compoundNames
makeMultiplePoolNames = function(compoundNames) {
  elementNames = mapply(function (cNames, eName) rep(eName, length(cNames)), compoundNames, names(compoundNames), SIMPLIFY = FALSE)
  PoolNames = mapply(makePoolName, compoundNames, elementNames, SIMPLIFY = F)
}

subsetGangstas = function(gangstaObjects, attributeName, attributeValue) {
  if(attributeName == "class") {
    itMatches = sapply(gangstaObjects, is, class2 = attributeValue)
  } else {
    if(is.logical(attributeValue)) {
      itMatches = sapply(gangstaObjects, "[", name=attributeName)
      ## where itMatches is TRUE, return TRUE; where NULL or FALSE, return FALSE
      itMatches = sapply(itMatches, function(x) !is.null(x) && x)
    } else {
      itMatches = (sapply(gangstaObjects, "[", name=attributeName) == attributeValue)
    }
  }
  return(gangstaObjects[itMatches])
}

gangstasExist = function(gangstaObjects, gangstaNames, checkClass = "") {
  ## The next line calls "getGangsta()", which will halt execution
  ## and throw and error if a gangsta doesn't exist.
  matchingGangstas = getGangstas(gangstaObjects, gangstaNames)
  if(any(checkClass!="")) {
    badClass = !sapply(matchingGangstas, is, checkClass)
    if(any(badClass)) {
      stop(paste0("Gangsta objects of type '", checkClass, "' are required for, but the following requested gangsta objects are not of that class: ", paste0(gangstaNames[badClass], collapse = ", ")))
    }
  }
  return(TRUE)
}

getGangstas = function(gangstaObjects, gangstaNames) {
  hits = lapply(gangstaNames, subsetGangstas, gangstaObjects = gangstaObjects, attributeName = "name")
  notFound = (sapply(hits, length) == 0)
  if(any(notFound)){
    stop(paste("Gangstas with the following names were requested but not found in the list of gangsta objects: ", paste0(gangstaNames[notFound], collapse = ", ")))
  }
  duplicates = (sapply(hits, length) > 1)
  if(any(duplicates)) {
    stop(paste0("More than one gangsta object exists with the following names ", paste0(gangstaNames[duplicates], collapse = ", ")))
  }
  return(unlist(hits, recursive = F))
}

getGangstaAttribute = function(gangstaObjects, attribName) {
  sapply(gangstaObjects, "[[", attribName)
}


# findByAttribute = function(gangstaList, attribName, value) {
#   itMatches = (sapply(gangstaList, "[", name=attribName) == value)
#   gangstaList[itMatches]
# }


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




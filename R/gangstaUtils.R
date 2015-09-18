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




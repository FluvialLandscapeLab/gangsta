isotopePostProcess <- function(results,
                               gangstas,
                               isotopesToTrack,
                               initialIsotopicRatios,
                               delValCalculators){
  # Find names of pools with isotope tracking
  elementNames = names(isotopesToTrack)
  poolObjects = subsetGangstas(gangstas, "class", getOption("gangsta.classes")["pool"])
  poolObjectsByElement = lapply(elementNames,
                                subsetGangstas,
                                gangstaObjects = poolObjects,
                                attributeName = "elementName")
  poolObjectsWithIsotopeTracking = unlist(poolObjectsByElement, recursive = F)
  poolNames = names(poolObjectsWithIsotopeTracking)
  poolObjectElements = getGangstaAttribute(poolObjectsWithIsotopeTracking, "elementName")

  # Get names of transfers into each pool with isotope tracking
  transfers = subsetGangstas(gangstas, "class", getOption("gangsta.classes")["trans"])
  transfersIn = lapply(poolNames,
                       subsetGangstas,
                       gangstaObjects = transfers,
                       attributeName = getOption("gangsta.attributes")["toPool"])
  transferInNames = lapply(transfersIn, getGangstaAttribute,
                           getOption("gangsta.attributes")["name"])
  transferInVars = lapply(transferInNames, makeTransferMolTransVars)
  names(transferInVars) = names(poolObjectsWithIsotopeTracking)

  # Get names of pools from which elements are transferred into each pool with isotope tracking
  transfersInFromPoolNames = lapply(transfersIn, getGangstaAttribute, getOption("gangsta.attributes")["fromPool"])

  # Get names of transfers out of each pool with isotope tracking
  transfersOut = lapply(names(poolObjectsWithIsotopeTracking),
                        subsetGangstas, gangstaObjects = transfers,
                        attributeName = getOption("gangsta.attributes")["fromPool"])
  transferOutNames = lapply(transfersOut, getGangstaAttribute,
                            getOption("gangsta.attributes")["name"])
  transferOutVars = lapply(transferOutNames, makeTransferMolTransVars)
  names(transferOutVars) = names(poolObjectsWithIsotopeTracking)

  # Order list of initial isotopic ratios according to pool objects by element list
  initialIsotopicRatios = initialIsotopicRatios[names(poolObjectsWithIsotopeTracking)]

  # Calculate atoms of isotope in each pool, atoms transferred in and out, and final isotopic composition
  # Calculate the number of atoms of isotope initially in each pool
  for(i in 1:length(results)){
    # Calculate the number of atoms of isotope initially in each pool
    if(i == 1){
      results[[i]]$isotopeVals$initialMols = mapply(
        function(poolNames, initialIsotopicRatios){
          initialIsotopicRatios * results[[1]]$poolVals[poolNames, "initial"]
        },
        poolNames = names(poolObjectsWithIsotopeTracking),
        initialIsotopicRatios = initialIsotopicRatios
      )
      results[[i]]$isotopeVals$initialIsotopicRatios = initialIsotopicRatios
    }else{
      results[[i]]$isotopeVals$initialMols =  results[[i-1]]$isotopeVals$finalMols
      results[[i]]$isotopeVals$initialIsotopicRatios = results[[i-1]]$isotopeVals$finalIsotopicRatios
    }

    # Get atoms transferred for transfers in and out of each pool
    transferAtomsIn = lapply(transferInVars,
                             function(j) {
                               if(length(j>=1)){results[[i]]$molTransferValsByProcess[j,"molTransfered"]}else{0}})
    names(transferAtomsIn) = names(poolObjectsWithIsotopeTracking)

    transferAtomsOut = lapply(transferOutVars,
                              function(x,i,j){
                                if(length(j>=1)){x[i,j]}else{0}},
                              i = 1, x = results)
    names(transferAtomsOut) = names(poolObjectsWithIsotopeTracking)


  }

  return(results)
}


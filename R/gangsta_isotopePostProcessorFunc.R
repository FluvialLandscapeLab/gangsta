isotopePostProcess = function(results,
                              gangstas,
                              isotopesToTrack,
                              initialIsotopicRatios,
                              isotopeLeakInList){
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
  transfersOut = lapply(poolNames,
                        subsetGangstas,
                        gangstaObjects = transfers,
                        attributeName = getOption("gangsta.attributes")["fromPool"])
  transferOutNames = lapply(transfersOut, getGangstaAttribute,
                            getOption("gangsta.attributes")["name"])
  transferOutVars = lapply(transferOutNames, makeTransferMolTransVars)
  names(transferOutVars) = names(poolObjectsWithIsotopeTracking)
  # Order list of initial isotopic ratios according to pool objects by element list
  initialIsotopicRatios = initialIsotopicRatios[names(poolObjectsWithIsotopeTracking)]
  # Calculate atoms of isotope in each pool, atoms transferred in and out, final atoms,  and final isotopic composition
  for(i in 1:length(results)){
    # Get mols of isotope leaked in
    stepIsotopeLeakInList = isotopeLeakInList[[i]]
    # Create list of zeros to initialize leak in list for all pools with isotope tracking
    poolMolsAddedByIsotope = lapply(poolObjectElements,
                                    function(element) rep(0, times = length(isotopesToTrack[[element]])))
    names(poolMolsAddedByIsotope) = poolNames

    # set the values of the initalized leak in list for any compounds that are non-zero
    if( length(stepIsotopeLeakInList [[1]]) > 0 ) {
      leakInPoolNames = mapply("[[", stepIsotopeLeakInList, 1)
      # poolMolsAddedByIsotope[leakInPoolNames] = unlist(mapply("[[", stepIsotopeLeakInList, 2, SIMPLIFY = FALSE) * results[[i]]$leakInPoolVals[leakInPoolNames])
      poolMolsAddedByIsotope[leakInPoolNames] = mapply(function(atomicFraction, leakInPoolName) {atomicFraction[[2]]*results[[i]]$leakInPoolVals[leakInPoolName]},
                                                       atomicFraction = stepIsotopeLeakInList,
                                                       leakInPoolName = leakInPoolNames,
                                                       SIMPLIFY = F )
    }
    # Calculate the number of atoms of isotope initially in each pool and add mols added
    if(i == 1){
      results[[i]]$isotopeVals$initialAtoms = mapply(
        function(poolNames, initialIsotopicRatios){
          initialIsotopicRatios * results[[1]]$poolVals[poolNames, "initial"]
        },
        poolNames = names(poolObjectsWithIsotopeTracking),
        initialIsotopicRatios = initialIsotopicRatios,
        SIMPLIFY = FALSE
      )
      results[[i]]$isotopeVals$initialIsotopicRatios = initialIsotopicRatios
    }else{
      results[[i]]$isotopeVals$initialAtoms =  mapply(function(existingAtomsIsotope, atomsIsotopeLeakedIn) existingAtomsIsotope + atomsIsotopeLeakedIn,
                                                      existingAtomsIsotope = results[[i-1]]$isotopeVals$finalAtoms,
                                                      atomsIsotopeLeakedIn = poolMolsAddedByIsotope,
                                                      SIMPLIFY = F)
      # results[[i]]$isotopeVals$initialIsotopicRatios = results[[i-1]]$isotopeVals$finalIsotopicRatios
      results[[i]]$isotopeVals$initialIsotopicRatios = mapply(function(isotopeAtoms, poolAtoms) if(poolAtoms>0){isotopeAtoms/poolAtoms}else{isotopeAtoms*NA},
                                                              isotopeAtoms = results[[i]]$isotopeVals$initialAtoms,
                                                              poolAtoms = results[[i]]$poolVals[poolNames,"initial"],
                                                              SIMPLIFY = F)
    }
    # Get atoms transferred for transfers in and out of each pool
    transferAtomsIn = lapply(transferInVars,
                             function(j) {
                               if(length(j>=1)){
                                 molsTransferred = results[[i]]$molTransferValsByProcess[j,"molTransfered"]
                                 noTransfer = is.na(molsTransferred)
                                 molsTransferred[noTransfer] = 0
                                 return(molsTransferred)
                               }else{0}
                             })
    transferAtomsOut = lapply(transferOutVars,
                              function(j) {
                                if(length(j>=1)){
                                  molsTransferred = results[[i]]$molTransferValsByProcess[j,"molTransfered"]
                                  noTransfer = is.na(molsTransferred)
                                  molsTransferred[noTransfer] = 0
                                  return(molsTransferred)
                                }else{0}
                              })
    # Get isotopic ratios for transfer pools
    transfersInIsotopicRatios = lapply(transfersInFromPoolNames,
                                       function(fromPoolNames){
                                         if(length(fromPoolNames)>=1){
                                           results[[i]]$isotopeVals$initialIsotopicRatios[fromPoolNames]
                                         }else{0}
                                       })
    transfersOutIsotopicRatios = lapply(poolNames,
                                        function(poolName)results[[i]]$isotopeVals$initialIsotopicRatios[poolName])
    names(transfersInIsotopicRatios) = poolNames
    names(transfersOutIsotopicRatios) = poolNames
    # Calculate net atoms transferred for each isotope in each pool according to
    # the transfer atoms and initial isotope ratios
    atomsTransferred = mapply(function(transferAtomsIn,
                                       transfersInIsotopicRatios,
                                       transferAtomsOut,
                                       transfersOutIsotopicRatios){
      atomsTransferredIn = mapply(function(isotopicRatios, transferAtoms) {
        isotopicRatios*transferAtoms
      },
      isotopicRatios = transfersInIsotopicRatios,
      transferAtoms = transferAtomsIn)

      atomsTransferredOut = mapply(function(isotopicRatios, transferAtoms){
        isotopicRatios*transferAtoms
      },
      transferAtoms = transferAtomsOut,
      isotopicRatios = transfersOutIsotopicRatios)

      if(!is.null(ncol(atomsTransferredIn))){
        totalAtomsTransferredInByIsotope = apply(X = atomsTransferredIn,
                                                 MARGIN = 1,
                                                 FUN = sum,
                                                 na.rm = T)
      }else{totalAtomsTransferredInByIsotope = 0}

      if(!is.null(ncol(atomsTransferredOut))){
        totalAtomsTransferredOutByIsotope = apply(X = atomsTransferredOut,
                                                  MARGIN = 1,
                                                  FUN = sum,
                                                  na.rm = T)
      }else{totalAtomsTransferredOutByIsotope = 0}
      return(totalAtomsTransferredInByIsotope - totalAtomsTransferredOutByIsotope)
    },
    transferAtomsIn = transferAtomsIn,
    transfersInIsotopicRatios = transfersInIsotopicRatios,
    transferAtomsOut = transferAtomsOut,
    transfersOutIsotopicRatios = transfersOutIsotopicRatios,
    SIMPLIFY = F
    )
    # Calculate final atoms of each isotope in each pool
    results[[i]]$isotopeVals$finalAtoms = mapply(function(initialAtoms, atomsTransferred){
      initialAtoms + atomsTransferred
    },
    initialAtoms = results[[i]]$isotopeVals$initialAtoms,
    atomsTransferred = atomsTransferred,
    SIMPLIFY = F)
    # Calculate isotopic ratios for each pool
    results[[i]]$isotopeVals$finalIsotopicRatios = mapply(function(isotopeAtoms, poolAtoms) if(poolAtoms>0){isotopeAtoms/poolAtoms}else{isotopeAtoms*NA},
                                                          isotopeAtoms = results[[i]]$isotopeVals$finalAtoms,
                                                          poolAtoms = results[[i]]$poolVals[poolNames,"final"],
                                                          SIMPLIFY = F)
  }
  return(results)
}

calculateDelVals = function(results, poolName, AtomicWeight, RStd){
  AtomicFractions = c(list(results[[1]]$isotopeVals$initialIsotopicRatios[[poolName]]),
                      lapply(results, function(timeStepResults){
                        timeStepResults$isotopeVals$finalIsotopicRatios[[poolName]]
                      }))
  names(AtomicFractions)[1]<-"Initial_Conditions"
  sapply(AtomicFractions, function(AtomicFractions){
    AF = AtomicFractions[AtomicWeight]
    R<- AF/(1-AF)
    del <- round((R/RStd-1)*1000, digits = 0)
  })
}

getPoolNamesForOneElement = function(results, gangstas, elementName){
  elementNames = names(isotopesToTrack)
  poolObjects = subsetGangstas(gangstas, "class", getOption("gangsta.classes")["pool"])
  poolObjectsOfElement = subsetGangstas(gangstaObjects = poolObjects,
                                        attributeName = "elementName",
                                        attributeValue = elementName)
  poolNames = names(poolObjectsOfElement)
  return(poolNames)
}

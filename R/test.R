gangstaCompounds = function(compoundParams) {
  return(unlist(lapply(compoundParams, do.call, what = "compoundFactory"), recursive = F))
}

gangstaProcesses = function(processParams, compounds) {

  #  expandedSpecs = unlist(lapply(processParams, expandMultiprocessSpec), recursive = F)
  #  processes = unlist(lapply(expandedSpecs, function(x) do.call(processFactory, args = c(list(compounds), x))), recursive = F)
  processes = unlist(lapply(processParams, function(x) do.call(processFactory, args = c(list(compounds), x))), recursive = F)
  names(processes) = getGangstaAttribute(processes, gangstaAttributeName("name"))
  return(processes)
}

gangstaBuildCompoundsAndProcesses = function(compoundParams, processParams){
  compounds = gangstaCompounds(compoundParams)
  processes = gangstaProcesses(processParams, compounds)
  return(c(compounds, processes))
}

leakIn = function(compoundNames){

  leakInList = list(
    list(
      list(compoundName = "Het", additionalMols = 1),
      list(compoundName = "Aut", additionalMols = 0),
      list(compoundName = "Met", additionalMols = 0),
      list(compoundName = "DOM", additionalMols = 0.1),
      list(compoundName = "O2" , additionalMols = 0.4)
    ),
    list(
      list()
    ),
    list(
      list(compoundName = "DOM", additionalMols = 0.1),
      list(compoundName = "NO3", additionalMols = 0.1),
      list(compoundName = "SO4" , additionalMols = 0.1)
    ),
    list(
      list(compoundName = "DOM", additionalMols = 0.3)
    ),
    list(
      list(compoundName = "O2" , additionalMols = 0.8)
    ),
    list(
      list(compoundName = "O2" , additionalMols = 0.8)
    ),
    list(
      list(compoundName = "DOM", additionalMols = 0.1),
      list(compoundName = "NO3", additionalMols = 0.1),
      list(compoundName = "SO4" , additionalMols = 0.1)
    ),
    list(
      list(compoundName = "DOM", additionalMols = 0.3)
    ),
    list(
      list(compoundName = "DOM", additionalMols = 0.3)
    )
  )

  leakInList = lapply(
    leakInList,
    function(x) {
      compoundExists = unlist(lapply(x, function(y) y$compoundName %in% compoundNames))
      if (!any(compoundExists)) {
        x = list(list())
      } else {
        x = x[compoundExists]
      }
      return(x)
    }
  )
  #   for (element in c("S", "O")) {
  #     if(!get(element)) {
  #       leakInList = lapply(
  #         leakInList,
  #         function(x) {
  #           if(is.null(names(x)))
  #           {
  #             return(x)
  #           } else {
  #             shortList = x[names(x)!=element]
  #             if(length(shortList)==0) {
  #               return(list(list()))
  #             } else {
  #               return(shortList)
  #             }
  #           }
  #         }
  #       )
  #     }
  #   }

  return(leakInList)
}

iterateGangsta = function(gangstaObjects, lpObject, leakInList = leakIn(), omitInitialRun = TRUE){
  startSuffix = gangstaVarName("startSuffix")
  gangstaResultsList = list()
  nameAttrName = gangstaAttributeName("name")
  numberOfSteps = length(leakInList) + 1

  gangstaCompounds = subsetGangstas(gangstaObjects,"class", "compound")
  gangstaCompoundNames = names(gangstaCompounds)
  #   gangstaPools = subsetGangstas(gangstaObjects,"class", "pool")
  #   gangstaPoolNames = names(gangstaPools)
  gangstaSinkCompounds = subsetGangstas(gangstaCompounds, "sourceSink", T)
  gangstaSinkCompoundNames = names(gangstaSinkCompounds)
  gangstaNonSinkCompounds = gangstaCompounds[!gangstaCompoundNames %in% gangstaSinkCompoundNames]
  gangstaNonSinkCompoundsNames = gangstaCompoundNames[!gangstaCompoundNames %in% gangstaSinkCompoundNames]
  #   gangstaSinkPoolIndex = getGangstaAttribute(gangstaPools, "compoundName") %in% gangstaSinkCompoundNames
  #   gangstaSinkPools = gangstaPools[gangstaSinkPoolIndex]
  #   gangstaSinkPoolNames = names(gangstaSinkPools)
  #   gangstaNonSinkPools = gangstaPools[!gangstaPoolNames %in% gangstaSinkPoolNames]
  #   gangstaNonSinkPoolsNames = gangstaPoolNames[!gangstaPoolNames %in% gangstaSinkPoolNames]

  listLoc = 0
  compoundMolsToAdd = NULL

  for(i in 1:numberOfSteps) {

    if(!omitInitialRun || i > 1) {
      lpStatus = suppressWarnings(lpSolveAPI::solve.lpExtPtr(lpObject))
      poolDifDF = poolDifs(gangstaObjects, lpObject, simple = FALSE)
      mt = massTransfers(gangstaObjects, lpObject, simple = FALSE, byProcess = FALSE)
      mtp = massTransfers(gangstaObjects, lpObject, simple = TRUE, byProcess = TRUE)
      pe = processEnergies(gangstaObjects, lpObject, simple = FALSE, catabolic = TRUE, anabolic = TRUE)
      compoundDifDF = compoundDifs(gangstaObjects, lpObject, simple = F)

      stepOutput = list(objective = lpSolveAPI::get.objective(lpObject),
                        status = lpStatus,
                        poolVals = poolDifDF,
                        massTransferVals = mt,
                        massTransferValsByProcess = mtp,
                        leakInVals = compoundMolsToAdd,
                        compoundVals = compoundDifDF,
                        processEnergyVals = pe)
      listLoc = listLoc + 1
      gangstaResultsList[[listLoc]] = stepOutput
      names(gangstaResultsList)[listLoc] = paste0("Iteration_", listLoc)

      if(i == numberOfSteps) break
    }

    stepLeakInList = unlist(leakInList[i], recursive = FALSE)
    names(stepLeakInList) = NULL
    numberOfLeakIns = length(stepLeakInList)

    if( length(stepLeakInList[[1]]) > 0 ) {
      compoundNames = mapply("[[", stepLeakInList, 1, SIMPLIFY = FALSE)
      compoundMolsToAdd = unlist(mapply("[[", stepLeakInList, 2, SIMPLIFY = FALSE))
      names(compoundMolsToAdd) = compoundNames

      compounds = gangstaCompounds[gangstaCompoundNames %in% compoundNames]
      #      pools = unlist(lapply(compoundNames, function(c) subsetGangstas(gangstaPools, "compoundName", c)), recursive = FALSE)
      #      poolNames = names(pools)

      #      refPoolNames = unlist(lapply(compoundNames, function(c) getGangstaAttribute(getGangstas(gangstaObjects, c), "referencePoolName")))
      #      names(refPoolMolsToAdd) = refPoolNames

      #      stoichPoolNames = poolNames[!(!is.na(match(poolNames, refPoolNames)))]
      #     } else {
      #       stoichPoolNames = NULL
      #       refPoolMolsToAdd = NULL
      #     }
      #     if(length(stoichPoolNames) > 0) {
      #       stoichMultipliers = getGangstaAttribute(pools[!is.na(match(poolNames, stoichPoolNames))], "molarRatio")
      #       stoichPools = pools[!is.na(match(poolNames, stoichPoolNames))]
      #       stoichPoolCompoundNames = getGangstaAttribute(stoichPools, "compoundName")
      #       stoichPoolCompounds = compounds[!is.na(match(compoundNames, stoichPoolCompoundNames))]
      #       stoichPoolReferencePoolNames = getGangstaAttribute(stoichPoolCompounds, "referencePoolName")
      #       refPoolMolsToMultiplyToStoichPoolMols = unlist(refPoolMolsToAdd[!is.na(match(refPoolNames, stoichPoolReferencePoolNames))])
      #       stoichPoolMolsToAdd = stoichMultipliers * refPoolMolsToMultiplyToStoichPoolMols
      #       molsToAdd = c(refPoolMolsToAdd, stoichPoolMolsToAdd)
      #     } else {
      #       molsToAdd = refPoolMolsToAdd
      #       if(length(molsToAdd) ==0) {
      #         molsToAdd = rep(0, length(gangstaNonSinkPoolsNames))
      #         names(molsToAdd) = gangstaNonSinkPoolsNames
      #       }
            if(length(compoundMolsToAdd) == 0) {
              compoundMolsToAdd = rep(0, length(gangstaNonSinkCompoundsNames))
              names(compoundMolsToAdd) = gangstaNonSinkCompoundsNames
            }
          }
      sinkCompoundWithLeakNames = names(compoundMolsToAdd)[names(compoundMolsToAdd) %in% gangstaSinkCompoundNames]
      if(length(sinkCompoundWithLeakNames) > 0) {
        stop("Elements are leaked into the following source/sinks: ", paste0(sinkCompoundWithLeakNames, collapse = "; "))
      }
      #     sinkPoolWithLeakNames = names(molsToAdd)[names(molsToAdd) %in% gangstaSinkPoolNames]
      #     if(length(sinkPoolWithLeakNames) > 0) {
      #       stop("Elements are leaked into the following source/sink pools: ", paste0(sinkPoolWithLeakNames, collapse = "; "))
      #     }

      if(!omitInitialRun || i>1) {
        # nextRunInitalVals = poolDifDF$final
        nextRunInitalVals = compoundDifDF$final
        names(nextRunInitalVals) = row.names(compoundDifDF)
        # names(nextRunInitalVals) = row.names(poolDifDF)
        nextRunInitalVals[gangstaSinkCompoundNames] = 0.0
        compoundMolsToAddIndices = match(names(compoundMolsToAdd), row.names(compoundDifDF))
        nextRunInitalVals[compoundMolsToAddIndices] = compoundMolsToAdd + nextRunInitalVals[compoundMolsToAddIndices]
      } else {
        nextRunInitalVals = compoundMolsToAdd
      }

      names(nextRunInitalVals) = paste(names(nextRunInitalVals), startSuffix, sep = ".")

      ### update initial mols using previous values + leak in list
      initCompoundIndexes = match(names(nextRunInitalVals), dimnames(lpObject)[[2]])
      lpSolveAPI::set.bounds(lpObject, lower = nextRunInitalVals, upper = nextRunInitalVals, columns = initCompoundIndexes)
      # initPoolIndexes = match(names(nextRunInitalVals), dimnames(lpObject)[[2]])
      # lpSolveAPI::set.bounds(lpObject, lower = nextRunInitalVals, upper = nextRunInitalVals, columns = initPoolIndexes)
    }
    return(gangstaResultsList)
  }

  doItGangsta = function(gangstaObjects, lpID, file = file.choose()){
    equations = makeEquations(gangstaObjects)
    writeGangstaModel(equations, file)
    gangsta.lp = readGangsta.lp(file)
    modelName = paste0("lp.", lpID)
    assign(modelName, gangsta.lp, envir = .GlobalEnv)
    return("Happy, happy!  Joy, joy!")
  }

  # doAll =  function(gangstaObjects = gangstas, file = "M:\\gangsta\\lpFiles\\test.lp"){
  #   modelName = doItGangsta(1, file)
  #   results = iterateGangsta(gangstaObjects, get(modelName, envir = .GlobalEnv))
  #   massTransfersPlot(gangstaObjects, results)
  # }

  doAll = function(tag, compoundParams, processParams, compoundNames) {
    gangstasName = paste0("gangstas", tag)
    resultsName = paste0("results", tag)
    modelName = paste0("lp.", tag)

    assign(gangstasName,
           gangstaBuildCompoundsAndProcesses(compoundParams, processParams),
           envir = .GlobalEnv)
    gangstaObjects = get(gangstasName, envir = .GlobalEnv)

    doItGangsta(gangstaObjects, tag, file = paste0("lpFiles/", tag, ".lp"))

    assign(resultsName,
           iterateGangsta(gangstaObjects, get(modelName, envir = .GlobalEnv), leakInList = leakIn(compoundNames)),
           envir = .GlobalEnv)
    massTransfersPlot(gangstaObjects, get(resultsName, envir = .GlobalEnv))
  }

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
  return(leakInList)
}

iterateGangsta = function(gangstaObjects, lpObject, leakInList = leakIn()){
  startSuffix = gangstaVarName("startSuffix")
  gangstaResultsList = list()
  nameAttrName = gangstaAttributeName("name")
  numberOfSteps = length(leakInList)

  gangstaCompounds = subsetGangstas(gangstaObjects,"class", "compound")
  gangstaCompoundNames = names(gangstaCompounds)
  gangstaSinkCompounds = subsetGangstas(gangstaCompounds, "sourceSink", T)
  gangstaSinkCompoundNames = names(gangstaSinkCompounds)

  initCompoundVarNames = makeCompoundStartMassVars(gangstaCompoundNames)
  initCompoundIndexes = match(initCompoundVarNames, dimnames(lpObject)[[2]])

  existingCompoundVals = lpSolveAPI::get.bounds(lpObject, columns = initCompoundIndexes)
  if(!(identical(existingCompoundVals$lower, existingCompoundVals$upper))) {
    stop("Upper and lower bounds for compound mols in existing lp models don't match.")
  }

  for(i in 1:numberOfSteps) {

    stepLeakInList = unlist(leakInList[i], recursive = FALSE)
    names(stepLeakInList) = NULL
    numberOfLeakIns = length(stepLeakInList)

    # create list of zeros to initials leak in list for all compounds
    compoundMolsAdded = rep(0, length(gangstaCompoundNames))
    names(compoundMolsAdded) = gangstaCompoundNames

    # set the values of the initalized leak in list for any compounds that are non-zero
    if( length(stepLeakInList[[1]]) > 0 ) {
      compoundNames = mapply("[[", stepLeakInList, 1)
      leakToSinks = compoundNames %in% gangstaSinkCompoundNames
      if(any(leakToSinks)) {
        stop("Compounds are leaked into the following source/sinks: ", paste0(compoundNames[leakToSinks], collapse = "; "))
      }
      compoundMolsAdded[compoundNames] = unlist(mapply("[[", stepLeakInList, 2, SIMPLIFY = FALSE))
    }

    ### update initial mols using leak in list
    if(i == 1) {
      existingCompoundVals = lpSolveAPI::get.bounds(lpObject, columns = initCompoundIndexes)$lower
    } else {
      existingCompoundVals =
        structure(
          gangstaResultsList[[i-1]]$compoundVals$final,
          names = row.names(gangstaResultsList[[i-1]]$compoundVals)
        )
    }

    initCompoundVals = compoundMolsAdded + existingCompoundVals

    lpSolveAPI::set.bounds(lpObject,
                           lower = initCompoundVals,
                           upper = initCompoundVals,
                           columns = initCompoundIndexes
    )

    lpStatus = suppressWarnings(lpSolveAPI::solve.lpExtPtr(lpObject))

    compoundDifDF = compoundDifs(gangstaObjects, lpObject, simple = F)
    poolDifDF = poolDifs(gangstaObjects, lpObject, simple = FALSE)
    if(i == 1) {
      ##  For the first timestep, I set the leak in values in the results to =
      ##  the sum of the initial values + the leak in list values
      poolMolsAdded = poolDifDF$initial
      compoundMolsAdded = initCompoundVals
    } else {
      poolMolsAdded = poolDifDF$initial - gangstaResultsList[[i-1]]$poolVals$final
      poolMolsAdded = ifelse(abs(poolMolsAdded) < 1e-12, 0, poolMolsAdded)
    }
    names(poolMolsAdded) = row.names(poolDifDF)

    mt = massTransfers(gangstaObjects, lpObject, simple = FALSE, byProcess = FALSE)
    mtp = massTransfers(gangstaObjects, lpObject, simple = TRUE, byProcess = TRUE)
    pe = processEnergies(gangstaObjects, lpObject, simple = FALSE, catabolic = TRUE, anabolic = TRUE)

    stepOutput = list(objective = lpSolveAPI::get.objective(lpObject),
                      status = lpStatus,
                      poolVals = poolDifDF,
                      massTransferVals = mt,
                      massTransferValsByProcess = mtp,
                      leakInPoolVals = poolMolsAdded,
                      leakInCompoundVals = compoundMolsAdded,
                      compoundVals = compoundDifDF,
                      processEnergyVals = pe)
    gangstaResultsList[[i]] = stepOutput
    names(gangstaResultsList)[i] = paste0("Iteration_", i)
  }
  return(gangstaResultsList)
}

doItGangsta = function(gangstaObjects, tag, file = file.choose()){
  equations = makeEquations(gangstaObjects)
  writeGangstaModel(equations, file)
  gangsta.lp = readGangsta.lp(file)
  modelName = paste0("lp.", tag)
  assign(modelName, gangsta.lp, envir = .GlobalEnv)
  return("Happy, happy!  Joy, joy!")
}

# doAll =  function(gangstaObjects = gangstas, file = "M:\\gangsta\\lpFiles\\test.lp"){
#   modelName = doItGangsta(1, file)
#   results = iterateGangsta(gangstaObjects, get(modelName, envir = .GlobalEnv))
#   massTransfersPlot(gangstaObjects, results)
# }

doAll = function(tag, compoundParams, processParams, compoundNames, sourceSinks) {
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
  massTransfersPlot(gangstaObjects, get(resultsName, envir = .GlobalEnv), tag = tag, sourceSinks = sourceSinks)
}

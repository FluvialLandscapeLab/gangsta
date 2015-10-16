
gangstaCompounds = function() {
  compoundParams = list(
    list(compoundName = "Het", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6),
    list(compoundName = "Aut", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6),
    list(compoundName = "Met", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6),
    list(compoundName = "DOM", molarRatios = c(C=1, N=16/106), initialMols = 0),
    list(compoundName = "XOM", molarRatios = c(C=1, N=16/106), initialMols = 0),
    list(compoundName = "CH4", molarRatios = c(C=1), initialMols = 0),
    list(compoundName = "NH4", molarRatios = c(N=1), initialMols = 0),
    list(compoundName = "NO3", molarRatios = c(N=1), initialMols = 0),
    #     list(compoundName = "NO2", molarRatios = c(N=1), initialMols = 0),
    #     list(compoundName = "N2O", molarRatios = c(N=1), initialMols = 0),
    list(compoundName = "O2" , molarRatios = c(O=1), initialMols = 0),
    list(compoundName = "SO4", molarRatios = c(S=1), initialMols = 0),
    list(compoundName = "CO2", molarRatios = c(C=1), initialMols = 0, sourceSink = T),
    list(compoundName = "N2" , molarRatios = c(N=1), initialMols = 0, sourceSink = T),
    list(compoundName = "HS" , molarRatios = c(S=1), initialMols = 0, sourceSink = T),
    list(compoundName = "Ox" , molarRatios = c(O=1), initialMols = 0, sourceSink = T)
  )
  compounds <<- unlist(lapply(compoundParams, do.call, what = "compoundFactory"), recursive = F)
}

gangstaProcesses = function() {
  processParams = list(
    list(
      name = "AssimDOM",
      energyTerm = -4.32E-04,
      fromCompoundNames = list(C = "DOM", N = "DOM"),
      toCompoundNames = list(C = ".", N = "."),
      molarTerms = list(C = 1, N = NA),
      organismName = "Het"
    ),
    list(
      name = "AssimCO2",
      energyTerm = -3.5E-02,
      fromCompoundNames = list(C = "CO2"),
      toCompoundNames = list(C = "."),
      molarTerms = list(C = 1),
      organismName = "Aut"
    ),
    list(
      name = "AssimCH4",
      energyTerm = -1.09E-03,
      fromCompoundNames = list(C = "CH4"),
      toCompoundNames = list(C = "."),
      molarTerms = list(C = 1),
      organismName = "Met"
    ),
    list(
      name = "AssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3"),
      toCompoundNames = list(N = "."),
      molarTerms = list(N = 1),
      organismName = c("Het", "Aut", "Met")
    ),
#     list(
#       name = "AssimNO2",
#       energyTerm = -1.25E-04,
#       fromCompoundNames = list(N = "NO2"),
#       toCompoundNames = list(N = "."),
#       molarTerms = list(N = 1),
#       organismName = c("Het", "Aut", "Met")
#     ),
    list(
      name = "AssimNH4",
      energyTerm = -3.18E-05,
      fromCompoundNames = list(N = "NH4"),
      toCompoundNames = list(N = "."),
      molarTerms = list(N = 1),
      organismName = c("Het", "Aut", "Met"),
      limitToInitMols = c(F, T, T)
    ),
    list(
      name = "Aerobic",
      energyTerm = 4.37-04,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), O = "O2"),
      toCompoundNames = list(C = "CO2", N = "NH4", O = "Ox"),
      molarTerms = list(C = 1, N = NA, O = 2),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "Denit",
      energyTerm = 2.88E-4 + 4.15E-4 + 6.45E-4,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "NO3"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2"),
      molarTerms = list(C = 3, N = NA, N = 8),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),

    list(
      name = "SulfateRed",
      energyTerm = 3.8E-05,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), S = "SO4"),
      toCompoundNames = list(C = "CO2", N = "NH4", S = "HS"),
      molarTerms = list(C = 1, N = NA, S = 0.5),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "Methanogenesis",
      energyTerm = 2.8E-05,
      fromCompoundNames = list(C = c(".", "DOM"), C = c(".", "DOM"), N = c(".", "DOM")),
      toCompoundNames = list(C = "CO2", C = "CH4", N = "NH4"),
      molarTerms = list(C = 0.5, C = 0.5, N = NA),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "Nitrif",
      energyTerm = 1.83E-04 + 1.48E-04,
      fromCompoundNames = list(N = "NH4", O = "O2"),
      toCompoundNames = list(N = "NO3", O = "Ox"),
      molarTerms = list(N = 8/3, O = 4),
      organismName = "Aut"
    ),
    list(
      name = "MethaneOxid",
      energyTerm = 4.09E-04,
      fromCompoundNames = list(C = "CH4", O = "O2"),
      toCompoundNames = list(C = "CO2", O = "Ox"),
      molarTerms = list(C = 0.5, O = 2),
      organismName = "Met"
    ),
    list(
      name = "Decay",
      energyTerm = 0,
      fromCompoundNames = list(C = ".", N = "."),
      toCompoundNames = list(C = "XOM", N = "XOM"),
      molarTerms = list(C = 1, N = NA),
      organismName = c("Aut", "Met")
    )
#     list(
#       name = "DenitStep1",
#       energyTerm = 2.88E-04,
#       fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "NO3"),
#       toCompoundNames = list(C = "CO2", N = "NH4", N = "NO2"),
#       molarTerms = list(C = 1, N = NA, N = 2),
#       organismName = "Het",
#       processSuffix = c("ofHet", "ofDOM")
#     ),
#     list(
#       name = "DenitStep2",
#       energyTerm = 4.15E-04,
#       fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "NO2"),
#       toCompoundNames = list(C = "CO2", N = "NH4", N = "N2O"),
#       molarTerms = list(C = 1, N = NA, N = 2),
#       organismName = "Het",
#       processSuffix = c("ofHet", "ofDOM")
#     ),
#     list(
#       name = "DenitStep3",
#       energyTerm = 6.45E-04,
#       fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "N2O"),
#       toCompoundNames = list(C = "CO2", N = "NH4", N = "N2"),
#       molarTerms = list(C = 1, N = NA, N = 4),
#       organismName = "Het",
#       processSuffix = c("ofHet", "ofDOM")
#     ),
#     list(
#       name = "NitrifStep1",
#       energyTerm = 1.83E-04,
#       fromCompoundNames = list(N = "NH4", O = "O2"),
#       toCompoundNames = list(N = "NO2", O = "Ox"),
#       molarTerms = list(N = 2/3, O = 2),
#       organismName = "Aut"
#     ),
#     list(
#       name = "NitrifStep2",
#       energyTerm = 1.48E-04,
#       fromCompoundNames = list(N = "NO2", O = "O2"),
#       toCompoundNames = list(N = "NO3", O = "Ox"),
#       molarTerms = list(N = 2, O = 2),
#       organismName = "Aut"
#     ),

  )

#   processParams = list(
#     list(
#       name = "Methanogenesis",
#       energyTerm = 2.8E-05,
#       fromCompoundNames = list(C = c(".", "DOM"), C = c(".", "."), N = c(".", "DOM"), N = c("NH4", ".")),
#       toCompoundNames = list(C = "CO2", C = "CH4", N = "NH4", N = c("NO3", "NH4")),
#       molarTerms = list(C = 0.5, C = 0.5, N = c(NA, NA), N = c(1, NA)),
#       organismName = "Het",
#       processSuffix = c("ofHet", "ofDOM")
#     )
#   )

  expandedSpecs = unlist(lapply(processParams, expandMultiprocessSpec), recursive = F)
  processes <<- unlist(lapply(expandedSpecs, function(x) do.call(processFactory, args = c(list(compounds), x))), recursive = F)
  names(processes) <<- getGangstaAttribute(processes, gangstaAttributeName("name"))
}

gangstaTest = function(){
  gangstaCompounds()
  gangstaProcesses()
  return(c(compounds, processes))
}

leakIn = function(){
  leakInList = list(
    list(
      list(compoundName = "Het", additionalMols = 1.1),
      list(compoundName = "Aut", additionalMols = 1.2),
      list(compoundName = "Met", additionalMols = 1.3),
      list(compoundName = "DOM", additionalMols = 1.4),
      list(compoundName = "O2" , additionalMols = 2)
    ),
    list(
      list(compoundName = "Aut", additionalMols = 1.5),
      list(compoundName = "NH4", additionalMols = 0.2),
      list(compoundName = "O2" , additionalMols = 1.4)
    ),
    list(
      list(compoundName = "O2" , additionalMols = 2)
    ),
    list(
      list(compoundName = "Met", additionalMols = 0.3),
      list(compoundName = "CH4", additionalMols = 0.6),
      list(compoundName = "O2" , additionalMols = 1.4)
    ),
    list(
      list(compoundName = "DOM", additionalMols = 2),
      list(compoundName = "NO3", additionalMols = 0.1),
      list(compoundName = "SO4", additionalMols = 1.0)
    ),
    list(
      list(compoundName = "Het", additionalMols = 0.1),
      list(compoundName = "Aut", additionalMols = 0.2),
      list(compoundName = "Met", additionalMols = 0.3),
      list(compoundName = "DOM", additionalMols = 0.4),
      list(compoundName = "XOM", additionalMols = 0.5),
      list(compoundName = "CH4", additionalMols = 0.6),
      list(compoundName = "NH4", additionalMols = 0.7),
      list(compoundName = "NO3", additionalMols = 0.8),
      #     list(compoundName = "NO2", additionalMols = 0),
      #     list(compoundName = "N2O", additionalMols = 0),
      list(compoundName = "O2" , additionalMols = 0.9),
      list(compoundName = "SO4", additionalMols = 1.0),
      list(compoundName = "CO2", additionalMols = 1.1),
      list(compoundName = "N2" , additionalMols = 1.2),
      list(compoundName = "HS" , additionalMols = 1.3)
      # list(compoundName = "Ox" , additionalMols = 1.4)
    )
  )
}

iterateGangsta = function(gangstaObjects, lpObject, leakInList = leakIn(), omitInitialRun = TRUE){
  startSuffix = gangstaVarName("startSuffix")
  gangstaResultsList = list()
  nameAttrName = gangstaAttributeName("name")
  numberOfSteps = length(leakInList) + 1

  gangstaCmpds = subsetGangstas(gangstas,"class", "compound")
  gangstaCmpdNames = names(gangstaCmpds)
  gangstaPools = subsetGangstas(gangstas,"class", "pool")
  gangstaPoolNames = names(gangstaPools)

  listLoc = 0
  molsToAdd = NULL

  for(i in 1:numberOfSteps) {

    if(!omitInitialRun || i > 1) {
      lpStatus = suppressWarnings(lpSolveAPI::solve.lpExtPtr(lpObject))
      pd = poolDifs(gangstaObjects, lpObject, simple = FALSE)
      mt = massTransfers(gangstaObjects, lpObject, simple = FALSE, byProcess = FALSE)

      stepOutput = list(objective = lpSolveAPI::get.objective(lpObject),
                        status = lpStatus,
                        poolVals = pd,
                        massTransferVals = mt,
                        leakInVals = molsToAdd)
      listLoc = listLoc + 1
      gangstaResultsList[[listLoc]] = stepOutput
      names(gangstaResultsList)[listLoc] = paste0("Iteration_", listLoc)

      if(i == numberOfSteps) break
    }

    stepLeakInList = unlist(leakInList[i], recursive = FALSE)
    numberOfLeakIns = length(stepLeakInList)
    cmpdNames = mapply("[[", stepLeakInList, 1, SIMPLIFY = FALSE)
    refPoolMolsToAdd = unlist(mapply("[[", stepLeakInList, 2, SIMPLIFY = FALSE))

    cmpds = gangstaCmpds[gangstaCmpdNames %in% cmpdNames]
    pools = unlist(lapply(cmpdNames, function(c) subsetGangstas(gangstaPools, "compoundName", c)), recursive = FALSE)
    poolNames = names(pools)

    refPoolNames = unlist(lapply(cmpdNames, function(c) getGangstaAttribute(getGangstas(gangstas, c), "referencePoolName")))
    names(refPoolMolsToAdd) = refPoolNames

    stoichPoolNames = poolNames[!(!is.na(match(poolNames, refPoolNames)))]
    if(length(stoichPoolNames) > 0) {
      stoichMultipliers = getGangstaAttribute(pools[!is.na(match(poolNames, stoichPoolNames))], "molarRatio")
      stoichPools = pools[!is.na(match(poolNames, stoichPoolNames))]
      stoichPoolCompoundNames = getGangstaAttribute(stoichPools, "compoundName")
      stoichPoolCompounds = cmpds[!is.na(match(cmpdNames, stoichPoolCompoundNames))]
      stoichPoolReferencePoolNames = getGangstaAttribute(stoichPoolCompounds, "referencePoolName")
      refPoolMolsToMultiplyToStoichPoolMols = unlist(refPoolMolsToAdd[!is.na(match(refPoolNames, stoichPoolReferencePoolNames))])
      stoichPoolMolsToAdd = stoichMultipliers * refPoolMolsToMultiplyToStoichPoolMols
      molsToAdd = c(refPoolMolsToAdd, stoichPoolMolsToAdd)
    } else {
      molsToAdd = refPoolMolsToAdd
    }

    if(!omitInitialRun || i>1) {
      molIndices = match(names(molsToAdd), row.names(pd))
      allMolsToAdd = molsToAdd + pd$final[molIndices]
    } else {
      allMolsToAdd = molsToAdd
    }

    names(allMolsToAdd) = paste(names(allMolsToAdd), startSuffix, sep = ".")

    ### update initial mols using previous values + leak in list
    initPoolIndexes = match(names(allMolsToAdd), dimnames(lpObject)[[2]])
    lpSolveAPI::set.bounds(lpObject, lower = allMolsToAdd, upper = allMolsToAdd, columns = initPoolIndexes)
  }
  return(gangstaResultsList)
}





doItGangsta = function(lpNum, file = file.choose()){
  gangstas <<- gangstaTest()
  equations <<- makeEquations(gangstas)
  writeGangstaModel(equations, file)
  gangsta.lp = readGangsta.lp(file)
  assign(paste0("test.lp.0", lpNum), gangsta.lp, envir = .GlobalEnv)
  # test.lp <<- readGangsta.lp(file)
  return(get(paste0("test.lp.0", lpNum)))
}

doAll =  function(gangstaObjects = gangstas, file = "M:\\gangsta\\lpFiles\\test.lp"){
  doItGangsta(1, file)
  results <<- iterateGangsta(gangstaObjects, test.lp.01)
  massTransfersPlot(gangstaObjects, results)
}


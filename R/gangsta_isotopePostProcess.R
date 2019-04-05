# Load Gangsta Rapper Functions
source("./R/gangsta_rapper_functions.R")

# Read lp Model
lpModel = read.lp("./vignettes/vignette.lp", verbose = "normal")

# 0. Run a time-stepping version of the model
time <- 1:3
solve(lpModel)

  # Make a matrix to store model results, where columns are variables and rows are time steps
modelVarNames = dimnames(lpModel)[[2]]
results <- matrix(ncol = length(get.variables(lpModel)),
                  nrow = length(time),
                  dimnames = list(time, modelVarNames))
results[1,] <- get.variables(lpModel)

  # Find names of variables with initial and final values
initValNames <- getValueNamesToUpdate(lpModel)$initialValueNames
finalValNames <- getValueNamesToUpdate(lpModel)$finalValueNames

for(i in 2:length(time)){
  newInitValues <- results[i-1,][finalValNames]
  # rename the final values with initial names
  names(newInitValues) = initValNames
  # set init values in lpmodel
  mapply(
    function(solvedValue, initValueColumnNumber){
      lpSolveAPI::set.bounds(
        lpModel,
        lower = solvedValue,
        upper = solvedValue,
        columns = initValueColumnNumber
      )
    },
    newInitValues,
    match(initValNames, dimnames(lpModel)[[2]])
  )
  # solve model for current time step
  lpSolveAPI::solve.lpExtPtr(lpModel)
  results[i,] = lpSolveAPI::get.variables(lpModel)
}

# 1. Find names of pools with isotope tracking
isotopesToTrack = list(C = c(12,13), O = c(16,17,18))
elementNames = names(isotopesToTrack)
poolObjects = subsetGangstas(myGangstas, "class", getOption("gangsta.classes")["pool"])
poolObjectsByElement = lapply(elementNames,
                              subsetGangstas,
                              gangstaObjects = poolObjects,
                              attributeName = "elementName")
poolObjectsWithIsotopeTracking = unlist(poolObjectsByElement, recursive = F)


# 2. Define initial isotopic ratios for each isotope in each pool
initialIsotopicRatios = list(CO2_C = c("12" = 0.9, "13" = 0.1),
                             DOM_C = c("12" = 0.9, "13" = 0.1),
                             CH4_C = c("12" = 0.9, "13" = 0.1),
                             Aut_C = c("12" = 0.9, "13" = 0.1),
                             Met_C = c("12" = 0.9, "13" = 0.1),
                             CO2_O = c("16" = 0.9, "17" = 0.05, "18" = 0.05),
                             NO3_O = c("16" = 0.9, "17" = 0.05, "18" = 0.05),
                             O2_O = c("16" = 0.9, "17" = 0.05, "18" = 0.05),
                             DOM_O = c("16" = 0.9, "17" = 0.05, "18" = 0.05),
                             Ox_O = c("16" = 0.9, "17" = 0.05, "18" = 0.05),
                             Aut_O = c("16" = 0.9, "17" = 0.05, "18" = 0.05),
                             Met_O = c("16" = 0.9, "17" = 0.05, "18" = 0.05)
)
  # order list of initial isotopic ratios according to pool objects by element list
initialIsotopicRatios = initialIsotopicRatios[names(poolObjectsWithIsotopeTracking)]

#   Make a matrix to store isotopic ratios
isotopicRatios = lapply(initialIsotopicRatios,
                        function(initialIsotopicRatios){
                          isotopicRatioMatrix = matrix(ncol = length(initialIsotopicRatios),
                                                       nrow = length(time),
                                                       dimnames = list(time, names(initialIsotopicRatios)))
                          isotopicRatioMatrix[1,] = initialIsotopicRatios
                          return(isotopicRatioMatrix)
                        })

# 3. Calculate initial atoms of each isotope in each pool
  poolInitialAtomsVarNames = makePoolStartMolVars(names(poolObjectsWithIsotopeTracking))
  initialAtomsIsotope = mapply(
    function(poolInitialAtomsVarNames, initialIsotopicRatios){
      initialIsotopicRatios * results[1,poolInitialAtomsVarNames]
    },
    poolInitialAtomsVarNames = poolInitialAtomsVarNames,
    initialIsotopicRatios = initialIsotopicRatios
  )
  names(initialAtomsIsotope) = names(poolObjectsWithIsotopeTracking)

# 4. Get names of transfers in and out of pools with isotope tracking
  transfers = subsetGangstas(myGangstas, "class", getOption("gangsta.classes")["trans"])
  transfersIn = lapply(names(poolObjectsWithIsotopeTracking),
                       subsetGangstas, gangstaObjects = transfers,
                       attributeName = getOption("gangsta.attributes")["toPool"])
  transferInNames = lapply(transfersIn, getGangstaAttribute,
                           getOption("gangsta.attributes")["name"])
  transferInVars = lapply(transferInNames, makeTransferMolTransVars)
  names(transferInVars) = names(poolObjectsWithIsotopeTracking)

  transfersOut = lapply(names(poolObjectsWithIsotopeTracking),
                        subsetGangstas, gangstaObjects = transfers,
                        attributeName = getOption("gangsta.attributes")["fromPool"])
  transferOutNames = lapply(transfersOut, getGangstaAttribute,
                            getOption("gangsta.attributes")["name"])
  transferOutVars = lapply(transferOutNames, makeTransferMolTransVars)
  names(transferOutVars) = names(poolObjectsWithIsotopeTracking)

# 5. Get atoms transferred for transfers identified in step 4
  transferAtomsIn = lapply(transferInVars,
                               function(x,i,j) {
                                 if(length(j>=1)){x[i,j]}else{0}},
                               i = 1, x = results)
  names(transferAtomsIn) = names(poolObjectsWithIsotopeTracking)
  transferAtomsOut = lapply(transferOutVars,
                                function(x,i,j){
                                  if(length(j>=1)){x[i,j]}else{0}},
                                i = 1, x = results)
  names(transferAtomsOut) = names(poolObjectsWithIsotopeTracking)

# 6. Get initial isotopic ratios corresponding to transfers in and out
transfersInFromPoolNames = lapply(transfersIn, getGangstaAttribute, getOption("gangsta.attributes")["fromPool"])
transfersInFromPoolInitialAtomsVarNames = lapply(transfersInFromPoolNames, makePoolStartMolVars)
transfersInFromPoolInitialAtoms = lapply(transfersInFromPoolInitialAtomsVarNames,
                                         function(initialAtomsVarNames){
                                           results[1,initialAtomsVarNames]
                                         })
transfersInIsotopicRatios = lapply(transfersInFromPoolNames,
                                           function(fromPoolNames){
                                             if(length(fromPoolNames)>=1){
                                               lapply(fromPoolNames, function(fromPoolName) isotopicRatios[[fromPoolName]][1,])
                                             }else{
                                               0
                                             }
                                           })
transfersOutIsotopicRatios = lapply(names(poolObjectsWithIsotopeTracking),
                                    function(poolName) isotopicRatios[[poolName]][1,])

# 7. Calculate net atoms transferred for each isotope in each pool according to
# the transfer atoms and initial isotope ratios
atomsTransferred = mapply(function(transferAtomsIn,
                                   transfersInIsotopicRatios,
                                   transferAtomsOut,
                                   transfersOutIsotopicRatios){
  atomsTransferredIn = mapply(function(isotopicRatios, transferAtoms) isotopicRatios*transferAtoms,
                               isotopicRatios = transfersInIsotopicRatios,
                               transferAtoms = transferAtomsIn)
  atomsTransferredOut = mapply(function(isotopicRatios, transferAtoms) isotopicRatios*transferAtoms,
                               transferAtoms = transferAtomsOut,
                               MoreArgs = list(isotopicRatios = transfersOutIsotopicRatios))
  if(!is.null(ncol(atomsTransferredIn))){
    totalAtomsTransferredInByIsotope = apply(X = atomsTransferredIn,
                                             MARGIN = 1,
                                             FUN = sum)
  }else{totalAtomsTransferredInByIsotope = 0}

  if(!is.null(ncol(atomsTransferredOut)>=1)){
    totalAtomsTransferredOutByIsotope = apply(X = atomsTransferredOut,
                                              MARGIN = 1,
                                              FUN = sum)
  }else{totalAtomsTransferredOutByIsotope = 0}
  return(totalAtomsTransferredInByIsotope - totalAtomsTransferredOutByIsotope)
},
transferAtomsIn = transferAtomsIn,
transfersInIsotopicRatios = transfersInIsotopicRatios,
transferAtomsOut = transferAtomsOut,
transfersOutIsotopicRatios = transfersOutIsotopicRatios
)

# 8. Calculate initial and final atoms of each isotope in each pool
finalAtoms = mapply(function(initialAtoms, atomsTransferred){
  initialAtoms + atomsTransferred
},
initialAtoms = initialAtomsIsotope,
atomsTransferred = atomsTransferred)

# 9. Calculate isotopic ratios for each pool
poolFinalAtomsVarNames = makePoolEndMolVars(names(poolObjectsWithIsotopeTracking))
finalIsotopicRatios = mapply(function(isotopeAtoms, poolAtoms) if(poolAtoms>0){isotopeAtoms/poolAtoms}else{NA},
       isotopeAtoms = finalAtoms,
       poolAtoms = results[1,poolFinalAtomsVarNames])

# 10. Repeat steps 5-9 for additional time steps






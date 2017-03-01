#### Make transfers table
vectorizeTransferInfo = function(gangstas, info){
  nTransfers = length(subsetGangstas(gangstas, "class", "transformation"))

  transferInfo =
    unlist(
      lapply(
        1:nTransfers,
        function(x) subsetGangstas(gangstas, "class", "transformation")[[x]][info]
      ),
      use.names = F
    )
  return(transferInfo)
}

makeTransferTable = function(gangstas) {

  specifyDecimal = function(x, y) format(round(x, y), nsmall=y)

  processNames = vectorizeTransferInfo(gangstas, "processName")
  # transferNames = vectorizeTransferInfo(gangstas, "name")
  fromNames = vectorizeTransferInfo(gangstas, "from")
  toNames = vectorizeTransferInfo(gangstas, "to")
  organismNames = c("Het","Aut","Met")
  for(i in 1:length(organismNames)){
    fromNames = gsub(organismNames[i], "Biomass", fromNames)
    toNames = gsub(organismNames[i], "Biomass", toNames)
  }

  molarTerms = vectorizeTransferInfo(gangstas, "molarTerm")
  roundedMolarTerms = specifyDecimal(molarTerms, 3)
  simpleMolarTerms = ifelse(roundedMolarTerms == specifyDecimal(1, 3), "", roundedMolarTerms)
  sepVec = ifelse(molarTerms == 1, "", " ")
  fromValues = paste0(simpleMolarTerms, sepVec, fromNames)
  toValues = paste0(simpleMolarTerms, sepVec, toNames)
  transferValues = paste(fromValues, "-->", toValues)
  jouleToMolRatios = scales::scientific(vectorizeTransferInfo(gangstas, "molarAffinity"), digits = 3)
  organismName = substr(processNames, 1, 3)
  processSimpleNames = substr(processNames, 4, nchar(processNames))
  element = substr(fromNames, nchar(fromNames), nchar(fromNames))

  processSorter = rep(0, length(processSimpleNames))
  processSorter =
    ifelse(
      processSimpleNames %in% processSimpleNames[grep("Assim", processSimpleNames)], 1,
      ifelse(
        processSimpleNames %in% processSimpleNames[grep("Decay", processSimpleNames)], 2, 3
      )
    )
  elementSorter = rep(0, length(element))
  elementSorter =
    ifelse(element == "C", 1,
           ifelse(element == "H", 2,
                  ifelse(element == "O", 3,
                         ifelse(element == "N", 4,
                                ifelse(element == "S", 5, NA
                                )
                         )
                  )
           )
    )
  organismSorter = rep(0, length(organismName))
  organismSorter =
    ifelse(organismName == "Met", 1,
           ifelse(organismName == "Aut", 2,
                  ifelse(organismName == "Het", 3, NA
                  )
           )
    )
  organismSorter =
    ifelse(processSorter == 3, organismSorter, 1)

  transferTable =
    data.frame(
      processSorter,
      elementSorter,
      organismSorter,
      organism = organismName,
      process = processSimpleNames,
      element = element,
      # transferNames,
      transfer = transferValues,
      umols = roundedMolarTerms,
      jouleToMolRatio = jouleToMolRatios
    )

  transferTable =
    plyr::ddply(transferTable, c("processSorter", "organismSorter", "process", "elementSorter", "transfer", "element", "umols", "jouleToMolRatio"), plyr::summarise,
                organisms = paste(unique(organism), collapse = ", ")
                )
  transferTable = transferTable[, -which(names(transferTable) %in% c("processSorter", "organismSorter", "elementSorter"))]

  return(transferTable)
}









getEnergies = function(results, energyType= "catabolic"){

  nResults = length(results)
  resultList = lapply(results, "[[", "processEnergyVals")

  for(i in 1:nResults){

    result = subset(resultList[[i]], resultList[[i]]$procType == "catabolic")
    result$process = row.names(result)
    result$iteration = i

    if(i == 1) {
      energyResults = result
    } else {
      energyResults = rbind(energyResults, result)
    }
  }
  row.names(energyResults) = 1:nrow(energyResults)
  modelID = deparse(substitute(results))
  energyResults$modelID = rep(modelID, nrow(energyResults))
  return(energyResults)
}

# result = get(results)[i]
# ls(pattern = "results")

summariseEnergies = function(energyResults) {
  summarizedResults = list(
    modelID = unique(energyResults$modelID),
    potentialProc = unique(energyResults$process),
    realizedProc = unique(
      energyResults$process[energyResults$energy > 0]
    )
  )
  return(summarizedResults)
}


getCatabEnergies = function(resultNames = ls(pattern = "results")){
  nScenarios = length(results)

  energyResultsList = list()
  for(i in 1:nScenarios) {
    result = get(resultNames[i])
    energies = getEnergies(result)
    summarizedEnergies = summariseEnergies(energies)
    energyResultsList = c(energyResultsList, summarizedEnergies)
  }
  return(energyResultsList)
  }


getScenarioEnergies = function(resultNames = ls(pattern = "results"),
                               results = get(ls(pattern = "results"))) {
  nResults = length(results)
  resultList = list()

  for(i in 1:nResults) {
    result = summariseEnergies(getEnergies(results[i]))
    resultList[i] = list(result$modelID[1],
                         c(result$potentialProc, result$realizedProc))
    }
}

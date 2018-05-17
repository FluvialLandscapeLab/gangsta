exprsnPaste = function(...) {
  if(any(unlist(lapply(list(...), function(x) length(x) == 0)))) {
    result = character(0)
  } else {
    result = paste(...)
  }
  return(result)
}

exprsnPaste0 = function(...) {
  return(exprsnPaste(..., sep = ""))
}

makeGenericVars = function(prefixes, varTag, separator = ".") {
  return(exprsnPaste0(prefixes, separator, gangstaVarName(varTag)))
}

makeProcessEnergyVars = function(processNames) {
  return(makeGenericVars(processNames, "energySuffixProcess"))
}

makeOrgansimEnergyVars = function(organismNames) {
  return(makeGenericVars(organismNames, "energySuffixOrganism"))
}

makeCompoundStartMolVars = function(compoundNames) {
  return(makeGenericVars(compoundNames, "startSuffixCompound"))
  # return(makeGenericVars(compoundNames, "startSuffix"))
}

makeCompoundEndMolVars = function(compoundNames) {
  return(makeGenericVars(compoundNames, "endSuffixCompound"))
  # return(makeGenericVars(compoundNames, "endSuffix"))
}

makePoolStartMolVars = function(poolNames) {
  return(makeGenericVars(poolNames, "startSuffixPool"))
  # return(makeGenericVars(poolNames, "startSuffix"))
}

makePoolEndMolVars = function(poolNames) {
  return(makeGenericVars(poolNames, "endSuffixPool"))
  # return(makeGenericVars(poolNames, "endSuffix"))
}

makeRespEnergyVars = function(organismNames) {
  return(makeGenericVars(organismNames, "respEnergy"))
}

makeTransferMolTransVars = function(transferNames) {
  return(makeGenericVars(transferNames,"transSuffix"))
}


makeLPSolveHeader = function(headerText, majorHeader = F) {
  if(majorHeader){
    header = paste("/********************", headerText, "********************/")
  } else {
    header = paste("/*****", headerText, "*****/")
  }
  return(header)
}

makeExpressions = function(gangstaObjects) {

  ## Get gangsta.option values


  endSuffixPool = gangstaVarName("endSuffixPool")
  endSuffixCompound = gangstaVarName("endSuffixCompound")

  energySuffixProcess = gangstaVarName("energySuffixProcess")
  energySuffixOrganism = gangstaVarName("energySuffixOrganism")

  startSuffixPool = gangstaVarName("startSuffixPool")
  startSuffixCompound = gangstaVarName("startSuffixCompound")

  transSuffix = gangstaVarName("transSuffix")

  respEnergyVarName = gangstaVarName("respEnergy")
  respRateVarName = gangstaVarName("respRate")

  organismClassName = gangstaClassName("org")
  processClassName = gangstaClassName("proc")
  transClassName = gangstaClassName("trans")
  poolClassName = gangstaClassName("pool")
  compoundClassName = gangstaClassName("comp")
  metabolicClassName = gangstaClassName("metab")

  nameAttrName = gangstaAttributeName("name")
  orgNameAttrName = gangstaAttributeName("orgName")
  fromPoolAttrName = gangstaAttributeName("fromPool")
  toPoolAttrName = gangstaAttributeName("toPool")
  respAttrName = gangstaAttributeName("respRate")
  procNameAttrName = gangstaAttributeName("procName")
  energyToMolsAttrName =  gangstaAttributeName("molarAffinity")
  limitToStartAttrName = gangstaAttributeName("limitToStartMols")
  molarRatioAttrName = gangstaAttributeName("molRatio")
  initialMolsAttrName = gangstaAttributeName("initialMolecules")
  compNameAttrName =  gangstaAttributeName("compName")
  InfiniteCompoundAttrName = gangstaAttributeName("infiniteCompound")
  energyTermAttrName = gangstaAttributeName("energy")
  elementAttrName = gangstaAttributeName("element")
  transOptionsAttrName = gangstaAttributeName("transOptions")

  exprsnMaxBiomass = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)

    organismEndMolVarNames = makeCompoundEndMolVars(organismNames)

    goalExprsnHeader =
      makeLPSolveHeader("Optimization function (Exprsn. 1)", F)
    totalBiomassHeader =
      makeLPSolveHeader("The Total_Biomass.finalMolecules equals the sum of the Organism.finalMolecules for each organism type (Exprsn. 2)", F)

    expressions =
      c(
        goalExprsnHeader,
        paste0("MAX: ", makeCompoundEndMolVars("Total_Biomass")),
        totalBiomassHeader,
        paste(makeCompoundEndMolVars("Total_Biomass"),"=", paste0(organismEndMolVarNames, collapse = " + "))
      )
    return(expressions)
  }

  exprsnBiosynthesis = function() {
    poolObjects = subsetGangstas(gangstaObjects, "class", poolClassName)

    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)

    organismPoolObjects =
      unlist(
        lapply(organismNames, subsetGangstas, gangstaObjects = poolObjects, attributeName = compNameAttrName)
      )

    organismPoolNames = getGangstaAttribute(organismPoolObjects, nameAttrName)

    transferObjects = subsetGangstas(gangstaObjects, "class", transClassName)

    transferInObjects =
      lapply(organismPoolNames, subsetGangstas, gangstaObjects = transferObjects, attributeName = toPoolAttrName)
    transferInNames = lapply(transferInObjects, getGangstaAttribute, nameAttrName)
    transferInVars = lapply(transferInNames, makeTransferMolTransVars)

  }

  exprsnAllowNegatives = function() {
    processes = subsetGangstas(gangstaObjects, "class", processClassName)
    processNames = getGangstaAttribute(processes, nameAttrName)
    energyVarNames = makeProcessEnergyVars(processNames)

    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)
    respEnergyVars = makeRespEnergyVars(organismNames)

    decayDissimAndAssimExpressions = c(exprsnPaste("-Inf <", c(energyVarNames), "< +Inf"))
    respExpressions = c(exprsnPaste("-Inf <", c(respEnergyVars), "< +Inf"))

    compounds = subsetGangstas(gangstaObjects, "class", compoundClassName)
    compounds = subsetGangstas(compounds, InfiniteCompoundAttrName, T)
    compoundNames = getGangstaAttribute(compounds, nameAttrName)
    compoundVarNames = makeCompoundEndMolVars(compoundNames)

    InfiniteCompoundCompoundExpressions = c(exprsnPaste("-Inf <", compoundVarNames, "< +Inf"))

    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    pools = lapply(compoundNames, subsetGangstas, gangstaObjects = pools, attributeName = compNameAttrName)
    pools = unlist(pools, recursive = F)
    poolNames = getGangstaAttribute(pools, nameAttrName)
    poolVarNames = c(makePoolEndMolVars(poolNames), makePoolStartMolVars(poolNames))

    InfiniteCompoundPoolExpressions = c(exprsnPaste("-Inf <", poolVarNames, "< +Inf"))

    InfiniteCompoundCompoundHeader =
      makeLPSolveHeader("For each InfiniteCompound, remove constraints on InfiniteCompound.finalMolecules", F)
    InfiniteCompoundPoolHeader =
      makeLPSolveHeader("For each Pool associated with an InfiniteCompound, remove constraints on Pool.finalAtoms", F)
    decayDissimAndAssimHeader =
      makeLPSolveHeader("For each Process, the Process.netEnergy is unconstrained", F)
    respHeader =
      makeLPSolveHeader("For each Organism type, the Organism.respirationEnergy is unconstrained", F)

    expressions = c(
      InfiniteCompoundCompoundHeader,
      InfiniteCompoundCompoundExpressions,

      InfiniteCompoundPoolHeader,
      InfiniteCompoundPoolExpressions,

      decayDissimAndAssimHeader,
      decayDissimAndAssimExpressions,

      respHeader,
      respExpressions
      )
    return(expressions)
  }

  exprsnInitialMols = function() {
    compounds = subsetGangstas(gangstaObjects, "class", compoundClassName)
    compoundNames = getGangstaAttribute(compounds, nameAttrName)
    initialCompoundMols = getGangstaAttribute(compounds, initialMolsAttrName)
    compoundStartMolVarNames = makeCompoundStartMolVars(compoundNames)

    initialMoleculesHeader =
      makeLPSolveHeader("Set FiniteCompound.initialMolecules & InfiniteCompound.initialMolecules", F)
    expressions =
      c(
        initialMoleculesHeader,
        paste(compoundStartMolVarNames, "=", initialCompoundMols)
      )
    return(expressions)
  }

  exprsnCompoundStoich = function() {
    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    poolNames = getGangstaAttribute(pools, nameAttrName)
    molarRatios = getGangstaAttribute(pools, molarRatioAttrName)

    # get compound name and then compound objects and then initial compound mols for each pool
    compoundNames = getGangstaAttribute(pools, compNameAttrName)
    compounds = getGangstas(gangstaObjects, compoundNames)

    # Starting stoicheometry expressions
    poolStartMolVarNames = makePoolStartMolVars(poolNames)
    compoundStartMolVarNames = makeCompoundStartMolVars(compoundNames)
    startStoichHeader =
      makeLPSolveHeader("For each elemental Pool, the Pool.initialAtoms must conform to Compound stoichiometry (Exprsn. 3)", F)
    expressions =
      c(
        startStoichHeader,
        paste(poolStartMolVarNames, "=", molarRatios, compoundStartMolVarNames)
        )

    # Ending stoicheometry expressions
    poolEndMolVarNames = makePoolEndMolVars(poolNames)
    compoundEndMolVarNames = makeCompoundEndMolVars(compoundNames)
    endStoichHeader =
      makeLPSolveHeader("For each elemental Pool, the Pool.finalAtoms must conform to Compound stoichiometry (Exprsn. 4)", F)
    expressions =
      c(
        expressions,
        endStoichHeader,
        paste(poolEndMolVarNames, "=", molarRatios, compoundEndMolVarNames)
      )

    return(expressions)
  }

  exprsnProcessEnergy = function() {

    metabolicProcesses = subsetGangstas(gangstaObjects, "class", metabolicClassName)
    organismNames = unique(getGangstaAttribute(metabolicProcesses, orgNameAttrName))

    processesByOrganism = lapply(organismNames, subsetGangstas, gangstaObjects = metabolicProcesses, attributeName = orgNameAttrName)
    processNamesByOrganism = lapply(processesByOrganism, getGangstaAttribute, attribName = nameAttrName)

    processEnergyVars = lapply(processNamesByOrganism, makeProcessEnergyVars)
    orgEnergyVars = makeOrgansimEnergyVars(organismNames)

    exprsnRightHandSides = sapply(processEnergyVars, paste, collapse = " + ")
    expressions = paste(orgEnergyVars, "=", exprsnRightHandSides)

    processEnergyHeader =
      makeLPSolveHeader("For each organism type, the total energy associated with assimilatory and dissimilatory Processes must equal the sum of the Process.netEnergy for all Processes (Exprsn. 9)", F)

    expressions = c(processEnergyHeader, expressions)

    return(expressions)
  }

  exprsnRespEnergy = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)
    biomassEndVarNames = makeCompoundEndMolVars(organismNames)

    respEnergyVarNames = makeRespEnergyVars(organismNames)
    respirationRates = getGangstaAttribute(organisms, respAttrName)

    expressions = paste(respEnergyVarNames, "=", respirationRates, biomassEndVarNames)

    respEnergyHeader =
      makeLPSolveHeader("Maintenance respiration energy is proportional to the biomass of each organism type (Exprsn. 8)", F)
    expressions =
      c(
        respEnergyHeader,
        expressions
      )
    return(expressions)
  }

  exprsnEnergyBalance = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)

    respEnergyVars = makeRespEnergyVars(organismNames)
    orgEnergyVars = makeOrgansimEnergyVars(organismNames)

    expressions = paste("0 =", respEnergyVars, "+", orgEnergyVars)

    energyBalHeader =
      makeLPSolveHeader("Organisms must meet their respiratory energy demands (Exprsn. 10)", F)

    expressions = c(energyBalHeader, expressions)

    return(expressions)
  }

  exprsnTransfer = function() {

    metabolics = subsetGangstas(gangstaObjects, "class", metabolicClassName)
    metabolicNames = getGangstaAttribute(metabolics, nameAttrName)
    transferOptions = getGangstaAttribute(metabolics, transOptionsAttrName)

    transfers = subsetGangstas(gangstaObjects, "class", transClassName)
    metabolicTransfers = lapply(metabolicNames, subsetGangstas, gangstaObjects = transfers, attributeName = procNameAttrName)

    metabolicTransferNames = lapply(metabolicTransfers, getGangstaAttribute, attribName = nameAttrName)
    metabolicTransferMolTransVars = lapply(metabolicTransferNames, makeTransferMolTransVars)

    metabolicProcessEnergyVars = makeProcessEnergyVars(metabolicNames)

    energyToMolsRatios = lapply(metabolicTransfers, getGangstaAttribute, attribName = energyToMolsAttrName)

    nestedExprsnTransfer = function(metabolicProcessEnergyVars, energyToMolsRatios, metabolicTransferMolTransVars, transferOptions) {
        return(
          sapply(
            transferOptions,
            function(idx) {
              paste(
                metabolicProcessEnergyVars,
                "=",
                paste(energyToMolsRatios[idx], metabolicTransferMolTransVars[idx], collapse = " + ")
              )
            }
          )
        )
    }

    expressions =
      unlist(
        mapply(
          nestedExprsnTransfer,
          metabolicProcessEnergyVars,
          energyToMolsRatios,
          metabolicTransferMolTransVars,
          transferOptions
        )
      )

    transferHeader =
      makeLPSolveHeader("The number of atoms transferred during each Process must be in accordance with the stoichiometry and chemical affinity of the Process (Exprsn. 5)", F)

    expressions =
      c(
        transferHeader,
        expressions
      )
    return(expressions)
  }

  exprsnMolBalance = function() {
    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    poolNames = getGangstaAttribute(pools, nameAttrName)
    poolStartMolVars = makePoolStartMolVars(poolNames)
    poolEndMolVars = makePoolEndMolVars(poolNames)

    transfers = subsetGangstas(gangstaObjects, "class", transClassName)

    transfersIn = lapply(poolNames, subsetGangstas, gangstaObjects = transfers, attributeName = toPoolAttrName)
    transferInNames = lapply(transfersIn, getGangstaAttribute, nameAttrName)
    transferInVars = lapply(transferInNames, makeTransferMolTransVars)

    transfersOut = lapply(poolNames, subsetGangstas, gangstaObjects = transfers, attributeName = fromPoolAttrName)
    transferOutNames = lapply(transfersOut, getGangstaAttribute, nameAttrName)
    transferOutVars = lapply(transferOutNames, makeTransferMolTransVars)

    transferInSumString = lapply(transferInVars, function(x) paste0(" + ", paste0(x, collapse = " + ")))
    transferOutSumString = lapply(transferOutVars, function(x) paste0(" - ", paste0(x, collapse = " - ")))

    transferInSumString[sapply(transfersIn, length) == 0] = ""
    transferOutSumString[sapply(transfersOut, length) == 0] = ""

    expressions = paste0(poolEndMolVars, " = ", poolStartMolVars, transferInSumString, transferOutSumString)

    molarBalHeader =
      makeLPSolveHeader("Elemental molar balance must be conserved (Exprsn. 6)", F)

    expressions =
      c(
        molarBalHeader,
        expressions
      )

    return(expressions)
  }

  exprsnLimitToStartingMol = function() {
    compounds = subsetGangstas(gangstaObjects, "class", compoundClassName)
    compounds = subsetGangstas(compounds, InfiniteCompoundAttrName, F)
    compoundNames = getGangstaAttribute(compounds, nameAttrName)

    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    pools = lapply(compoundNames, subsetGangstas, gangstaObjects = pools, attributeName = compNameAttrName)
    pools = unlist(pools, recursive = F)

    poolNames = getGangstaAttribute(pools, nameAttrName)
    poolStartMolVars = makePoolStartMolVars(poolNames)

    transfers = subsetGangstas(gangstaObjects, "class", transClassName)

    transfersOut = lapply(poolNames, subsetGangstas, gangstaObjects = transfers, attributeName = fromPoolAttrName)
    # remove any transfers that user defines as not limited by the starting mol of the from pool.
    transfersOut = lapply(transfersOut, subsetGangstas, attributeName = limitToStartAttrName, attributeValue = T)

    poolStartMolVars = poolStartMolVars[sapply(transfersOut, length) > 0]
    transfersOut = transfersOut[sapply(transfersOut, length) > 0]

    transferOutNames = lapply(transfersOut, getGangstaAttribute, nameAttrName)
    transferOutVars = lapply(transferOutNames, makeTransferMolTransVars)

    transferOutSumString = lapply(transferOutVars, paste0, collapse = " + ")

    expressions = paste0(poolStartMolVars, " >= ", transferOutSumString)

    limToInitHeader =
      makeLPSolveHeader("Generally, Transfers out of each Pool must be less than the Pool.initialAtoms (Exprsn. 7)", F)
    expressions = c(limToInitHeader, expressions)

    return(expressions)
  }

  # goal Expressions
  goalExprsns = exprsnMaxBiomass()

  # Allow negative values where appropriate
  allowNegativeExprsns = exprsnAllowNegatives()

  # Initial Mols
  initialMolExprsns = exprsnInitialMols()

  # Compound Stoichiometry Expressions
  compoundStoichExprsns = exprsnCompoundStoich()

  # Respiration Expressions
  respEnergyExprsns = exprsnRespEnergy()

  # Process Energy Expressions
  processEnergyExprsns = exprsnProcessEnergy()

  # Energy Balance Expressions
  energyBalExprsns = exprsnEnergyBalance()

  # Transfer Expressions
  transferExprsns = exprsnTransfer()

  # MolBalance Expressions
  molBalExprsns = exprsnMolBalance()

  limitToStartExprsns = exprsnLimitToStartingMol()

  allExpressions = c(
    makeLPSolveHeader("OPTIMIZATION", T),
    goalExprsns,

    makeLPSolveHeader("INITIALIZATION", T),
    initialMolExprsns,
    allowNegativeExprsns,

    makeLPSolveHeader("STOICHIOMETRY", T),
    compoundStoichExprsns,

    makeLPSolveHeader("MOLAR TRANSFER", T),
    transferExprsns,

    makeLPSolveHeader("ELEMENTAL MOLAR BALANCE", T),
    molBalExprsns,
    limitToStartExprsns,

    makeLPSolveHeader("ENERGETIC CONSTRAINTS", T),
    respEnergyExprsns,
    processEnergyExprsns,
    energyBalExprsns
  )
  return(allExpressions)
}

formatExpressions = function(expressions) {
  expressions =
    sapply(
      expressions,
      function(eq)
        ifelse(
          substr(eq, 1, 2) != "/*",
          paste0(eq, ";"),
          eq
        )
    )
  return(expressions)
}

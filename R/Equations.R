eqnPaste = function(...) {
  if(any(unlist(lapply(list(...), function(x) length(x) == 0)))) {
    result = character(0)
  } else {
    result = paste(...)
  }
  return(result)
}

eqnPaste0 = function(...) {
  return(eqnPaste(..., sep = ""))
}

makeGenericVars = function(prefixes, varTag, separator = ".") {
  return(eqnPaste0(prefixes, separator, gangstaVarName(varTag)))
}

makeProcessEnergyVars = function(processNames) {
  return(makeGenericVars(processNames, "energySuffixProcess"))
}

makeOrgansimEnergyVars = function(organismNames) {
  return(makeGenericVars(organismNames, "energySuffixOrganism"))
}

makeCompoundStartMassVars = function(compoundNames) {
  return(makeGenericVars(compoundNames, "startSuffixCompound"))
  # return(makeGenericVars(compoundNames, "startSuffix"))
}

makeCompoundEndMassVars = function(compoundNames) {
  return(makeGenericVars(compoundNames, "endSuffixCompound"))
  # return(makeGenericVars(compoundNames, "endSuffix"))
}

makePoolStartMassVars = function(poolNames) {
  return(makeGenericVars(poolNames, "startSuffixPool"))
  # return(makeGenericVars(poolNames, "startSuffix"))
}

makePoolEndMassVars = function(poolNames) {
  return(makeGenericVars(poolNames, "endSuffixPool"))
  # return(makeGenericVars(poolNames, "endSuffix"))
}

makeRespEnergyVars = function(organismNames) {
  return(makeGenericVars(organismNames, "respEnergy"))
}

makeTransformationMassTransVars = function(transformationNames) {
  return(makeGenericVars(transformationNames,"transSuffix"))
}


makeLPSolveHeader = function(headerText, majorHeader = F) {
  if(majorHeader){
    header = paste("/********************", headerText, "********************/")
  } else {
    header = paste("/*****", headerText, "*****/")
  }
  return(header)
}

makeEquations = function(gangstaObjects) {

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
  # finalMoleculesAttrName = gangstaAttributeName("finalMolecules")
  compNameAttrName =  gangstaAttributeName("compName")
  InfiniteCompoundAttrName = gangstaAttributeName("InfiniteCompound")
  energyTermAttrName = gangstaAttributeName("energy")
  elementAttrName = gangstaAttributeName("element")
  transOptionsAttrName = gangstaAttributeName("transOptions")

  eqnMaxBiomass = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)

    organismEndMassVarNames = makeCompoundEndMassVars(organismNames)

    goalEqnHeader =
      makeLPSolveHeader("Optimization function (Exprsn. 1)", F)
    totalBiomassHeader =
      makeLPSolveHeader("The Total_Biomass.finalMolecules equals the sum of the Organism.finalMolecules for each organism type (Exprsn. 2)", F)

    equations =
      c(
        goalEqnHeader,
        paste0("MAX: ", makeCompoundEndMassVars("Total_Biomass")),
        totalBiomassHeader,
        paste(makeCompoundEndMassVars("Total_Biomass"),"=", paste0(organismEndMassVarNames, collapse = " + "))
      )
    return(equations)
  }

  eqnBiosynthesis = function() {
    poolObjects = subsetGangstas(gangstaObjects, "class", poolClassName)

    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)

    organismPoolObjects =
      unlist(
        lapply(organismNames, subsetGangstas, gangstaObjects = poolObjects, attributeName = compNameAttrName)
      )

    organismPoolNames = getGangstaAttribute(organismPoolObjects, nameAttrName)

    transformationObjects = subsetGangstas(gangstaObjects, "class", transClassName)

    transferInObjects =
      lapply(organismPoolNames, subsetGangstas, gangstaObjects = transformationObjects, attributeName = toPoolAttrName)
    transferInNames = lapply(transferInObjects, getGangstaAttribute, nameAttrName)
    transferInVars = lapply(transferInNames, makeTransformationMassTransVars)

  }

  eqnAllowNegatives = function() {
    processes = subsetGangstas(gangstaObjects, "class", processClassName)
    # negativeEnergyTerms = getGangstaAttribute(processes, energyTermAttrName) < 0
    # processes = processes[negativeEnergyTerms]
    processNames = getGangstaAttribute(processes, nameAttrName)
    energyVarNames = makeProcessEnergyVars(processNames)

    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)
    respEnergyVars = makeRespEnergyVars(organismNames)

    # respAndAssimEquations = c(eqnPaste("-Inf <", c(respEnergyVars, energyVarNames), "<= 0"))
    decayDissimAndAssimEquations = c(eqnPaste("-Inf <", c(energyVarNames), "< +Inf"))
    respEquations = c(eqnPaste("-Inf <", c(respEnergyVars), "< +Inf"))


    compounds = subsetGangstas(gangstaObjects, "class", compoundClassName)
    compounds = subsetGangstas(compounds, InfiniteCompoundAttrName, T)
    compoundNames = getGangstaAttribute(compounds, nameAttrName)
    compoundVarNames = makeCompoundEndMassVars(compoundNames)

    InfiniteCompoundCompoundEquations = c(eqnPaste("-Inf <", compoundVarNames, "< +Inf"))

    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    pools = lapply(compoundNames, subsetGangstas, gangstaObjects = pools, attributeName = compNameAttrName)
    pools = unlist(pools, recursive = F)
    poolNames = getGangstaAttribute(pools, nameAttrName)
    poolVarNames = c(makePoolEndMassVars(poolNames), makePoolStartMassVars(poolNames))

    InfiniteCompoundPoolEquations = c(eqnPaste("-Inf <", poolVarNames, "< +Inf"))

    InfiniteCompoundCompoundHeader =
      makeLPSolveHeader("For each InfiniteCompound, remove constraints on InfiniteCompound.finalMolecules", F)
    InfiniteCompoundPoolHeader =
      makeLPSolveHeader("For each Pool associated with an InfiniteCompound, remove constraints on Pool.finalAtoms", F)
    decayDissimAndAssimHeader =
      makeLPSolveHeader("For each Process, the Process.netEnergy is unconstrained", F)
    respHeader =
      makeLPSolveHeader("For each Organism type, the Organism.respirationEnergy is unconstrained", F)

    equations = c(
      InfiniteCompoundCompoundHeader,
      InfiniteCompoundCompoundEquations,

      InfiniteCompoundPoolHeader,
      InfiniteCompoundPoolEquations,

      decayDissimAndAssimHeader,
      decayDissimAndAssimEquations,

      respHeader,
      respEquations
      )
    return(equations)
  }

  eqnInitialMasses = function() {
    compounds = subsetGangstas(gangstaObjects, "class", compoundClassName)
    compoundNames = getGangstaAttribute(compounds, nameAttrName)
    initialCompoundMols = getGangstaAttribute(compounds, initialMolsAttrName)
    compoundStartMassVarNames = makeCompoundStartMassVars(compoundNames)

    initialMoleculesHeader =
      makeLPSolveHeader("Set FiniteCompound.initialMolecules & InfiniteCompound.initialMolecules", F)
    equations =
      c(
        initialMoleculesHeader,
        paste(compoundStartMassVarNames, "=", initialCompoundMols)
      )
    return(equations)
  }

  eqnCompoundStoich = function() {
    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    poolNames = getGangstaAttribute(pools, nameAttrName)
    molarRatios = getGangstaAttribute(pools, molarRatioAttrName)

    # get compound name and then compound objects and then initial compound mols for each pool
    compoundNames = getGangstaAttribute(pools, compNameAttrName)
    compounds = getGangstas(gangstaObjects, compoundNames)

    # Starting stoicheometry equations
    poolStartMassVarNames = makePoolStartMassVars(poolNames)
    compoundStartMassVarNames = makeCompoundStartMassVars(compoundNames)
    startStoichHeader =
      makeLPSolveHeader("For each elemental Pool, the Pool.initialAtoms must conform to Compound stoichiometry (Exprsn. 3)", F)
    equations =
      c(
        startStoichHeader,
        paste(poolStartMassVarNames, "=", molarRatios, compoundStartMassVarNames)
        )

    # Ending stoicheometry equations
    poolEndMassVarNames = makePoolEndMassVars(poolNames)
    compoundEndMassVarNames = makeCompoundEndMassVars(compoundNames)
    endStoichHeader =
      makeLPSolveHeader("For each elemental Pool, the Pool.finalAtoms must conform to Compound stoichiometry (Exprsn. 4)", F)
    equations =
      c(
        equations,
        endStoichHeader,
        paste(poolEndMassVarNames, "=", molarRatios, compoundEndMassVarNames)
      )

    return(equations)
  }

  eqnProcessEnergy = function() {

    metabolicProcesses = subsetGangstas(gangstaObjects, "class", metabolicClassName)
    organismNames = unique(getGangstaAttribute(metabolicProcesses, orgNameAttrName))

    processesByOrganism = lapply(organismNames, subsetGangstas, gangstaObjects = metabolicProcesses, attributeName = orgNameAttrName)
    processNamesByOrganism = lapply(processesByOrganism, getGangstaAttribute, attribName = nameAttrName)

    processEnergyVars = lapply(processNamesByOrganism, makeProcessEnergyVars)
    orgEnergyVars = makeOrgansimEnergyVars(organismNames)

    eqnRightHandSides = sapply(processEnergyVars, paste, collapse = " + ")
    equations = paste(orgEnergyVars, "=", eqnRightHandSides)

    processEnergyHeader =
      makeLPSolveHeader("For each organism type, the total energy associated with assimilatory and dissimilatory Processes must equal the sum of the Process.netEnergy for all Processes (Exprsn. 9)", F)

    equations = c(processEnergyHeader, equations)

    return(equations)
  }

  eqnRespEnergy = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)
    biomassEndVarNames = makeCompoundEndMassVars(organismNames)

    respEnergyVarNames = makeRespEnergyVars(organismNames)
    respirationRates = getGangstaAttribute(organisms, respAttrName)

    equations = paste(respEnergyVarNames, "=", respirationRates, biomassEndVarNames)

    respEnergyHeader =
      makeLPSolveHeader("Maintenance respiration energy is proportional to the biomass of each organism type (Exprsn. 8)", F)
    equations =
      c(
        respEnergyHeader,
        equations
      )
    return(equations)
  }

  eqnEnergyBalance = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)

    respEnergyVars = makeRespEnergyVars(organismNames)
    orgEnergyVars = makeOrgansimEnergyVars(organismNames)

    equations = paste("0 =", respEnergyVars, "+", orgEnergyVars)

    energyBalHeader =
      makeLPSolveHeader("Organisms must meet their respiratory energy demands (Exprsn. 10)", F)

    equations = c(energyBalHeader, equations)

    return(equations)
  }

  eqnTransformation = function() {

    metabolics = subsetGangstas(gangstaObjects, "class", metabolicClassName)
    metabolicNames = getGangstaAttribute(metabolics, nameAttrName)
    transferOptions = getGangstaAttribute(metabolics, transOptionsAttrName)

    transformations = subsetGangstas(gangstaObjects, "class", transClassName)
    metabolicTransformations = lapply(metabolicNames, subsetGangstas, gangstaObjects = transformations, attributeName = procNameAttrName)
  #  metabolicTransformations = unlist(metabolicTransformations, recursive = F)

    metabolicTransformationNames = lapply(metabolicTransformations, getGangstaAttribute, attribName = nameAttrName)
    metabolicTransformationMassTransVars = lapply(metabolicTransformationNames, makeTransformationMassTransVars)

  #  metabolicProcessNames = lapply(metabolicTransformations, getGangstaAttribute, attribName = procNameAttrName)
    metabolicProcessEnergyVars = makeProcessEnergyVars(metabolicNames)

    energyToMolsRatios = lapply(metabolicTransformations, getGangstaAttribute, attribName = energyToMolsAttrName)

    nestedEqnTransformation = function(metabolicProcessEnergyVars, energyToMolsRatios, metabolicTransformationMassTransVars, transferOptions) {
        return(
          sapply(
            transferOptions,
            function(idx) {
              paste(
                metabolicProcessEnergyVars,
                "=",
                paste(energyToMolsRatios[idx], metabolicTransformationMassTransVars[idx], collapse = " + ")
              )
            }
          )
        )
    }

    equations =
      unlist(
        mapply(
          nestedEqnTransformation,
          metabolicProcessEnergyVars,
          energyToMolsRatios,
          metabolicTransformationMassTransVars,
          transferOptions
        )
      )

    # equations = paste(metabolicProcessEnergyVars, "=", energyToMolsRatios, metabolicTransformationMassTransVars)

    transferHeader =
      makeLPSolveHeader("The number of atoms transferred during each Process must be in accordance with the stoichiometry and chemical affinity of the Process (Exprsn. 5)", F)

    equations =
      c(
        transferHeader,
        equations
      )
    return(equations)
  }

  eqnMassBalance = function() {
    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    poolNames = getGangstaAttribute(pools, nameAttrName)
    poolStartMassVars = makePoolStartMassVars(poolNames)
    poolEndMassVars = makePoolEndMassVars(poolNames)

    transformations = subsetGangstas(gangstaObjects, "class", transClassName)

    transfersIn = lapply(poolNames, subsetGangstas, gangstaObjects = transformations, attributeName = toPoolAttrName)
    transferInNames = lapply(transfersIn, getGangstaAttribute, nameAttrName)
    transferInVars = lapply(transferInNames, makeTransformationMassTransVars)

    transfersOut = lapply(poolNames, subsetGangstas, gangstaObjects = transformations, attributeName = fromPoolAttrName)
    transferOutNames = lapply(transfersOut, getGangstaAttribute, nameAttrName)
    transferOutVars = lapply(transferOutNames, makeTransformationMassTransVars)

    transferInSumString = lapply(transferInVars, function(x) paste0(" + ", paste0(x, collapse = " + ")))
    transferOutSumString = lapply(transferOutVars, function(x) paste0(" - ", paste0(x, collapse = " - ")))

    transferInSumString[sapply(transfersIn, length) == 0] = ""
    transferOutSumString[sapply(transfersOut, length) == 0] = ""

    equations = paste0(poolEndMassVars, " = ", poolStartMassVars, transferInSumString, transferOutSumString)

    molarBalHeader =
      makeLPSolveHeader("Elemental molar balance must be conserved (Exprsn. 6)", F)

    equations =
      c(
        molarBalHeader,
        equations
      )

    return(equations)
  }

  eqnLimitToStartingMass = function() {
    compounds = subsetGangstas(gangstaObjects, "class", compoundClassName)
    compounds = subsetGangstas(compounds, InfiniteCompoundAttrName, F)
    compoundNames = getGangstaAttribute(compounds, nameAttrName)

    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    pools = lapply(compoundNames, subsetGangstas, gangstaObjects = pools, attributeName = compNameAttrName)
    pools = unlist(pools, recursive = F)

    poolNames = getGangstaAttribute(pools, nameAttrName)
    poolStartMassVars = makePoolStartMassVars(poolNames)

    transformations = subsetGangstas(gangstaObjects, "class", transClassName)

    transfersOut = lapply(poolNames, subsetGangstas, gangstaObjects = transformations, attributeName = fromPoolAttrName)
    # remove any transformations that user defines as not limited by the starting mass of the from pool.
    transfersOut = lapply(transfersOut, subsetGangstas, attributeName = limitToStartAttrName, attributeValue = T)

    poolStartMassVars = poolStartMassVars[sapply(transfersOut, length) > 0]
    transfersOut = transfersOut[sapply(transfersOut, length) > 0]

    transferOutNames = lapply(transfersOut, getGangstaAttribute, nameAttrName)
    transferOutVars = lapply(transferOutNames, makeTransformationMassTransVars)

    transferOutSumString = lapply(transferOutVars, paste0, collapse = " + ")

    equations = paste0(poolStartMassVars, " >= ", transferOutSumString)

    limToInitHeader =
      makeLPSolveHeader("Generally, Transfers out of each Pool must be less than the Pool.initialAtoms (Exprsn. 7)", F)
    equations = c(limToInitHeader, equations)

    return(equations)
  }

  # goal Equations
  goalEqns = eqnMaxBiomass()

  # Allow negative values where appropriate
  allowNegativeEqns = eqnAllowNegatives()

  # Initial Masses
  initialMassEqns = eqnInitialMasses()

  # Compound Stoichiometry Equations
  compoundStoichEqns = eqnCompoundStoich()

  # Respiration Equations
  respEnergyEqns = eqnRespEnergy()

  # Process Energy Equations
  processEnergyEqns = eqnProcessEnergy()

  # Energy Balance Equations
  energyBalEqns = eqnEnergyBalance()

  # Transformation Equations
  transformationEqns = eqnTransformation()

  # MassBalance Equations
  massBalEqns = eqnMassBalance()

  limitToStartEqns = eqnLimitToStartingMass()

  allEquations = c(
    makeLPSolveHeader("OPTIMIZATION", T),
    goalEqns,

    makeLPSolveHeader("INITIALIZATION", T),
    initialMassEqns,
    allowNegativeEqns,

    makeLPSolveHeader("STOICHIOMETRY", T),
    compoundStoichEqns,

    makeLPSolveHeader("MOLAR TRANSFER", T),
    transformationEqns,

    makeLPSolveHeader("ELEMENTAL MOLAR BALANCE", T),
    massBalEqns,
    limitToStartEqns,

    makeLPSolveHeader("ENERGETIC CONSTRAINTS", T),
    respEnergyEqns,
    processEnergyEqns,
    energyBalEqns
  )
  return(allEquations)
}

formatEquations = function(equations) {
  equations =
    sapply(
      equations,
      function(eq)
        ifelse(
          substr(eq, 1, 2) != "/*",
          paste0(eq, ";"),
          eq
        )
    )
  return(equations)
}

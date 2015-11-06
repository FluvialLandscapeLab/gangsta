makeProcessEnergyVars = function(processNames) {
  return(paste0(processNames, ".", gangstaVarName("energySuffix")))
}

makeOrgansimEnergyVars = function(organismNames) {
  return(paste0(organismNames, ".", gangstaVarName("energySuffix")))
}

makePoolStartMassVars = function(poolNames) {
  return(paste0(poolNames, ".", gangstaVarName("startSuffix")))
}

makePoolEndMassVars = function(poolNames) {
  return(paste0(poolNames, ".", gangstaVarName("endSuffix")))
}

makeRespEnergyVars = function(organismNames) {
  return(paste0(organismNames, ".", gangstaVarName("respEnergy")))
}

makeTransformationMassTransVars = function(transformationNames) {
  return(paste0(transformationNames, ".", gangstaVarName("transSuffix")))
}

makeBiomassRemainingVars = function(poolNames) {
  return(paste0(poolNames, ".", gangstaVarName("endSuffix")))
}

makeEquations = function(gangstaObjects) {

  ## Get gangsta.option values
  endSuffix = gangstaVarName("endSuffix")
  energySuffix = gangstaVarName("energySuffix")
  massSuffix = gangstaVarName("massSuffix")
  startSuffix = gangstaVarName("startSuffix")
  transSuffix = gangstaVarName("transSuffix")

  respEnergyVarName = gangstaVarName("respEnergy")
  respRateVarName = gangstaVarName("respRate")

  organismClassName = gangstaClassName("org")
  processClassName = gangstaClassName("proc")
  transClassName = gangstaClassName("trans")
  poolClassName = gangstaClassName("pool")
  boundClassName = gangstaClassName("bnd")
  compoundClassName = gangstaClassName("comp")
  metabolicClassName = gangstaClassName("metab")

  nameAttrName = gangstaAttributeName("name")
  orgNameAttrName = gangstaAttributeName("orgName")
  refPoolAttrName = gangstaAttributeName("refPool")
  fromPoolAttrName = gangstaAttributeName("fromPool")
  toPoolAttrName = gangstaAttributeName("toPool")
  respAttrName = gangstaAttributeName("respRate")
  procNameAttrName = gangstaAttributeName("procName")
  energyToMolsAttrName =  gangstaAttributeName("joulesToMols")
  limitToStartAttrName = gangstaAttributeName("limitToStartMass")
  molarRatioAttrName = gangstaAttributeName("molRatio")
  initialMolsAttrName = gangstaAttributeName("initMols")
  compNameAttrName =  gangstaAttributeName("compName")
  sourceSinkAttrName = gangstaAttributeName("sourceSink")
  energyTermAttrName = gangstaAttributeName("energy")
  elementAttrName = gangstaAttributeName("element")

  eqnMaxBiomass = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)
    refPoolNames = getGangstaAttribute(organisms, refPoolAttrName)
    refPools = getGangstas(gangstaObjects, refPoolNames)
    refPoolElementNames = getGangstaAttribute(refPools, elementAttrName)
    if(length(unique(refPoolElementNames)) != 1) {
      stop("Reference pool elements must be the same across organisms.  Currently, reference pool elements are:\n", paste0("    ", organismNames, ": ", refPoolElementNames, collapse = "\n"))
    }
    refPoolEndMassVars = makePoolEndMassVars(refPoolNames)
    equations = c("MAX: Total_Biomass", paste0("Total_Biomass = ", paste0(refPoolEndMassVars, collapse = " + ")))
    return(equations)
  }

  eqnAllowNegatives = function() {
    processes = subsetGangstas(gangstaObjects, "class", processClassName)
    negativeEnergyTerms = getGangstaAttribute(processes, energyTermAttrName) < 0
    processes = processes[negativeEnergyTerms]
    processNames = getGangstaAttribute(processes, nameAttrName)
    energyVarNames = makeProcessEnergyVars(processNames)

    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)
    respEnergyVars = makeRespEnergyVars(organismNames)

    equations = paste("-Inf <", c(respEnergyVars, energyVarNames), "<= 0")

    compounds = subsetGangstas(gangstaObjects, "class", compoundClassName)
    compounds = subsetGangstas(compounds, sourceSinkAttrName, T)
    compoundNames = getGangstaAttribute(compounds, nameAttrName)

    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    pools = lapply(compoundNames, subsetGangstas, gangstaObjects = pools, attributeName = compNameAttrName)
    pools = unlist(pools, recursive = F)
    poolNames = getGangstaAttribute(pools, nameAttrName)
    poolVarNames = makePoolEndMassVars(poolNames)

    equations = c(equations, paste("-Inf <", poolVarNames, "< +Inf"))
    return(equations)
  }

  eqnInitialMasses = function() {
    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
    poolNames = getGangstaAttribute(pools, nameAttrName)
    compoundNames = getGangstaAttribute(pools, compNameAttrName)
    compounds = getGangstas(gangstaObjects, compoundNames)
    initialMols = getGangstaAttribute(compounds, initialMolsAttrName)
    isBound = sapply(pools, is, class2 = boundClassName)

    poolStartMassVars = makePoolStartMassVars(poolNames)

    equations = paste(poolStartMassVars[!isBound], "=", initialMols[!isBound])

    if(any(isBound)) {
      molarRatios = getGangstaAttribute(pools[isBound], molarRatioAttrName)
      refPoolNames = getGangstaAttribute(compounds[isBound], refPoolAttrName)
      refPoolStartMassVars = makePoolStartMassVars(refPoolNames)
      equations = c(equations, paste(poolStartMassVars[isBound], "=", molarRatios, refPoolStartMassVars))

      boundPoolEndMassVars = makePoolEndMassVars(poolNames[isBound])
      refPoolEndMassVars = makePoolEndMassVars(refPoolNames)
      equations = c(equations, paste(boundPoolEndMassVars, "=", molarRatios, refPoolEndMassVars))
    }

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
    return(equations)
  }

  eqnRespEnergy = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)
    referencePoolNames = getGangstaAttribute(organisms, refPoolAttrName)

    respEnergyVars = makeRespEnergyVars(organismNames)
    respirationRates = getGangstaAttribute(organisms, respAttrName)
    biomassVars = makeBiomassRemainingVars(referencePoolNames)

    equations = paste(respEnergyVars, "=", respirationRates, biomassVars)
    return(equations)
  }

  eqnEnergyBalance = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)

    respEnergyVars = makeRespEnergyVars(organismNames)
    orgEnergyVars = makeOrgansimEnergyVars(organismNames)
    return(paste("0 =", respEnergyVars, "+", orgEnergyVars))
  }

  eqnTransformation = function() {

    metabolics = subsetGangstas(gangstaObjects, "class", metabolicClassName)
    metabolicNames = getGangstaAttribute(metabolics, nameAttrName)

    transformations = subsetGangstas(gangstaObjects, "class", transClassName)
    metabolicTransformations = lapply(metabolicNames, subsetGangstas, gangstaObjects = transformations, attributeName = procNameAttrName)
    metabolicTransformations = unlist(metabolicTransformations, recursive = F)

    metabolicTransformationNames = getGangstaAttribute(metabolicTransformations, nameAttrName)
    metabolicTransformationMassTransVars = makeTransformationMassTransVars(metabolicTransformationNames)

    metabolicProcessNames = getGangstaAttribute(metabolicTransformations, procNameAttrName)
    metabolicProcessEnergyVars = makeProcessEnergyVars(metabolicProcessNames)

    energyToMolsRatios = getGangstaAttribute(metabolicTransformations, energyToMolsAttrName)

    equations = paste(metabolicProcessEnergyVars, "=", energyToMolsRatios, metabolicTransformationMassTransVars)
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
    return(equations)
  }

  eqnLimitToStartingMass = function() {
    compounds = subsetGangstas(gangstaObjects, "class", compoundClassName)
    compounds = subsetGangstas(compounds, sourceSinkAttrName, F)
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

    return(equations)
  }

  # goal Equations
  goalEqns = eqnMaxBiomass()

  # Allow negative values where appropriate
  allowNegativeEqns = eqnAllowNegatives()

  # Initial Masses
  initialMassEqns = eqnInitialMasses()

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
    goalEqns,
    initialMassEqns,
    allowNegativeEqns,
    respEnergyEqns,
    processEnergyEqns,
    energyBalEqns,
    transformationEqns,
    massBalEqns,
    limitToStartEqns
  )

  print("Damn it feels good to be a GANGSTA!")

  return(allEquations)
}


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

  nameAttrName = gangstaAttributeName("name")
  orgNameAttrName = gangstaAttributeName("orgName")
  refPoolAttrName = gangstaAttributeName("refPool")
  fromPoolAttrName = gangstaAttributeName("fromPool")
  toPoolAttrName = gangstaAttributeName("toPool")
  respAttrName = gangstaAttributeName("respRate")
  procNameAttrName = gangstaAttributeName("procName")
  energyToMassAttrName =  gangstaAttributeName("energyToMass")
  limitToStartAttrName = gangstaAttributeName("limitToStartMass")
  molarRatioAttrName = gangstaAttributeName("molRatio")
  initialMolsAttrName = gangstaAttributeName("initMols")
  compNameAttrName =  gangstaAttributeName("compName")

  makeProcessEnergyVars = function(processNames) {
    return(paste0(processNames, ".", energySuffix))
  }

  makeOrgansimEnergyVars = function(organismNames) {
    return(paste0(organismNames, ".", energySuffix))
  }

  makePoolStartMassVars = function(poolNames) {
    return(paste0(poolNames, ".", startSuffix))
  }

  makePoolEndMassVars = function(poolNames) {
    return(paste0(poolNames, ".", endSuffix))
  }

  makeRespEnergyVars = function(organismNames) {
    return(paste0(organismNames, ".", respEnergyVarName))
  }

  makeTransformationMassTransVars = function(transformationNames) {
    return(paste0(transformationNames, ".", transSuffix))
  }

  makeBiomassRemainingVars = function(poolNames) {
    return(paste0(poolNames, ".", endSuffix))
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
      refPoolNames = getGangstaAttribute(compounds[isBound], nameAttrName)
      refPoolStartMassVars = makePoolStartMassVars(refPoolNames)
      equations = c(equations, paste(poolStartMassVars[isBound], "=", molarRatios, refPoolStartMassVars))

      boundPoolEndMassVars = makePoolEndMassVars(poolNames[isBound])
      refPoolEndMassVars = makePoolEndMassVars(refPoolNames)
      equations = c(equations, paste(boundPoolEndMassVars, "=", molarRatios, refPoolEndMassVars))
    }

    return(equations)
  }

  eqnProcessEnergy = function() {

    processes = subsetGangstas(gangstaObjects, "class", processClassName)
    organismNames = unique(getGangstaAttribute(processes, orgNameAttrName))

    processesByOrganism = lapply(organismNames, subsetGangstas, gangstaObjects = processes, attributeName = orgNameAttrName)
    processNamesByOrganism = lapply(processesByOrganism, getGangstaAttribute, attribName = nameAttrName)

    processEnergyVars = lapply(processNamesByOrganism, makeProcessEnergyVars)
    orgEnergyVars <<- makeOrgansimEnergyVars(organismNames)

    eqnRightHandSides = sapply(processEnergyVars, paste, collapse = " + ")
    equations = paste(orgEnergyVars, "=", eqnRightHandSides)
    return(equations)
  }

  eqnRespEnergy = function() {
    organisms = subsetGangstas(gangstaObjects, "class", organismClassName)
    organismNames = getGangstaAttribute(organisms, nameAttrName)
    referencePoolNames = getGangstaAttribute(organisms, refPoolAttrName)

    respEnergyVars <<- makeRespEnergyVars(organismNames)
    respirationRates = getGangstaAttribute(organisms, respAttrName)
    biomassVars = makeBiomassRemainingVars(referencePoolNames)

    equations = paste(respEnergyVars, "=", respirationRates, biomassVars)
    return(equations)
  }

  eqnEnergyBalance = function() {
    return(paste("0 =", respEnergyVars, "+", orgEnergyVars))
  }

  eqnTransformation = function() {

    transformations = subsetGangstas(gangstaObjects, "class", transClassName)
    transformationNames = getGangstaAttribute(transformations, nameAttrName)
    transformationMassVars = makeTransformationMassTransVars(transformationNames)

    processNames = getGangstaAttribute(transformations, procNameAttrName)
    processEnergyVars = makeProcessEnergyVars(processNames)

    energyToMassRatios = getGangstaAttribute(transformations, energyToMassAttrName)

    equations = paste(transformationMassVars, "=", energyToMassRatios, processEnergyVars)
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
    pools = subsetGangstas(gangstaObjects, "class", poolClassName)
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
    initialMassEqns,
    respEnergyEqns,
    processEnergyEqns,
    energyBalEqns,
    transformationEqns,
    massBalEqns,
    limitToStartEqns
  )

  return(allEquations)
}


## ----eval = T------------------------------------------------------------
library(gangsta)

## ----eval = T------------------------------------------------------------
EcoStoichC = 106/106
EcoStoichN = 16/106
EcoStoichO = 110/106

## ----eval = T------------------------------------------------------------
resp = -2.83E-6

## ----eval = T------------------------------------------------------------
timeStep = 24

## ----eval = T------------------------------------------------------------
myGangstas =
  c(
    compoundFactory(
      compoundName = "CO2",
      molarRatios = c(C=1, O=2),
      initialMolecules = 0.2
    ),
    compoundFactory(
      compoundName = "CH4" ,
      molarRatios = c(C=1),
      initialMolecules = 0.1
      ),
    compoundFactory(
      compoundName = "NH4",
      molarRatios = c(N=1),
      initialMolecules = 0.3
    ),
    compoundFactory(
      compoundName = "NO3",
      molarRatios = c(N=1, O=3),
      initialMolecules = 0.1
      ),
    compoundFactory(
      compoundName = "O2",
      molarRatios = c(O=2),
      initialMolecules = 0.9
      ),
    compoundFactory(
      compoundName = "DOM",
      molarRatios = c(C=EcoStoichC, N=EcoStoichN, O=EcoStoichO),
      initialMolecules = 0.3
    )
  )

## ----eval = T------------------------------------------------------------
myGangstas =
  c(myGangstas,
    compoundFactory(
      compoundName = "Ox",
      molarRatios = c(O=1),
      initialMolecules = 0,
      infiniteCompound = T
      )
  )

## ----eval = T------------------------------------------------------------
myGangstas =
  c(myGangstas,
    compoundFactory(
      compoundName = "Aut",
      molarRatios = c(C=EcoStoichC, N=EcoStoichN, O=EcoStoichO),
      initialMolecules = 1,
      respirationRate = resp * timeStep
    ),
    compoundFactory(
      compoundName = "Met",
      molarRatios = c(C=EcoStoichC, N=EcoStoichN, O=EcoStoichO),
      initialMolecules = 1,
      respirationRate = resp * timeStep
    )
  )

## ----eval = T------------------------------------------------------------
myGangstas =
  c(
    myGangstas,
    processFactory(
      myGangstas,
      processName = "AutNitrif",
      energyTerm = 3.485E-4,
      fromCompoundNames = list(N = "NH4", O = "O2", O = "O2"),
      toCompoundNames = list(N = "NO3", O = "NO3", O = "Ox"),
      molarTerms = list(N = 1, O = 3, O = 1),
      organismName = "Aut"
    ),
    processFactory(
      myGangstas,
      processName = "MetMethaneOxid",
      energyTerm = 8.18E-4,
      fromCompoundNames = list(C = "CH4", O = "O2", O = "O2"),
      toCompoundNames = list(C = "CO2", O = "CO2", O = "Ox"),
      molarTerms = list(C = 1, O = 2, O = 2),
      organismName = "Met"
    )
  )

## ----eval = T------------------------------------------------------------
myGangstas =
  c(
    myGangstas,
    processFactory(
      myGangstas,
      processName = "AutAssimCO2",
      energyTerm = -3.5E-03,
      fromCompoundNames = list(C = "CO2", O = "CO2", O = "CO2"),
      toCompoundNames = list(C = "Aut", O = "Aut", O = "Ox"),
      molarTerms = list(C = 1, O = 2, O = 2),
      transferOptions = list(C = 1, O = 2:3),
      organismName = "Aut"
    ),
    processFactory(
      myGangstas,
      processName = "MetAssimCH4",
      energyTerm = -1.09E-03,
      fromCompoundNames = list(C = "CH4"),
      toCompoundNames = list(C = "Met"),
      molarTerms = list(C = 1),
      organismName = "Met"
    ),
    processFactory(
      myGangstas,
      processName = "AutAssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3", O = "NO3", O = "NO3"),
      toCompoundNames = list(N = "Aut", O = "Aut", O = "Ox"),
      molarTerms = list(N = 1, O = 3, O = 3),
      transferOptions = list(N = 1, O = 2:3),
      organismName = "Aut"
    ),
    processFactory(
      myGangstas,
      processName = "MetAssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3", O = "NO3", O = "NO3"),
      toCompoundNames = list(N = "Met", O = "Met", O= "Ox"),
      molarTerms = list(N = 1, O = 3, O = 3),
      transferOptions = list(N = 1, O = 2:3),
      organismName = "Met"
    ),
    processFactory(
      myGangstas,
      processName = "AutAssimNH4",
      energyTerm = -3.18E-05,
      fromCompoundNames = list(N = "NH4"),
      toCompoundNames = list(N = "Aut"),
      molarTerms = list(N = 1),
      organismName = "Aut"
    ),
    processFactory(
      myGangstas,
      processName = "MetAssimNH4",
      energyTerm = -3.18E-05,
      fromCompoundNames = list(N = "NH4"),
      toCompoundNames = list(N = "Met"),
      molarTerms = list(N = 1),
      organismName = "Met"
    )
  )

## ----eval = T------------------------------------------------------------
myGangstas =
  c(myGangstas,
    processFactory(
      myGangstas,
      processName = "AutDecay",
      energyTerm = 0,
      fromCompoundNames = list(C = "Aut", N = "Aut", O = "Aut"),
      toCompoundNames = list(C = "DOM", N = "DOM", O = "DOM"),
      molarTerms = list(C = 1, N = EcoStoichN, O = EcoStoichO),
      organismName = c("Aut")
    ),
    processFactory(
      myGangstas,
      processName = "MetDecay",
      energyTerm = 0,
      fromCompoundNames = list(C = "Met", N = "Met", O = "Met"),
      toCompoundNames = list(C = "DOM", N = "DOM", O = "DOM"),
      molarTerms = list(C = 1, N = EcoStoichN, O = EcoStoichO),
      organismName = c("Met")
    )
  )

## ----eval = FALSE--------------------------------------------------------
writeGangstaModel(gangstaObjects = myGangstas, file = file.choose())

## ----eval = FALSE--------------------------------------------------------
library(lpSolveAPI)

## ----eval = FALSE--------------------------------------------------------
lpModel = read.lp(file.choose(), verbose = "normal")

## ----eval = FALSE--------------------------------------------------------
solve(lpModel)

## ----eval = FALSE--------------------------------------------------------
results = get.variables(lpModel)
names(results) = dimnames(lpModel)[[2]]
results


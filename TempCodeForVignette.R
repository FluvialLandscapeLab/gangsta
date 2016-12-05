BioStoichN = 16/106
BioStoichO = 110/106
DOMStoichN = 16/106
DOMStoichO = 110/106


myGangstas =
  c(
    compoundFactory(
      compoundName = "Aut",
      molarRatios = c(C=1, N=BioStoichN, O=BioStoichO),
      initialMols = 1,
      respirationRate = -2.83E-6 * 24
    ),
    compoundFactory(
      compoundName = "Met",
      molarRatios = c(C=1, N=BioStoichN, O=BioStoichO),
      initialMols = 1,
      respirationRate = -2.83E-6 * 24
    ),
    compoundFactory(
      compoundName = "DOM",
      molarRatios = c(C=1, N=DOMStoichN, O=DOMStoichO),
      initialMols = 0.3
      ),
    compoundFactory(
      compoundName = "CH4" ,
      molarRatios = c(C=1),
      initialMols = 0.1
      ),
    compoundFactory(
      compoundName = "NH4",
      molarRatios = c(N=1),
      initialMols = 0.3
    ),
    compoundFactory(
      compoundName = "NO3",
      molarRatios = c(N=1, O=3),
      initialMols = 0.1
      ),
    compoundFactory(
      compoundName = "O2",
      molarRatios = c(O=2),
      initialMols = 0.9
      ),
    compoundFactory(
      compoundName = "CO2",
      molarRatios = c(C=1, O=2),
      initialMols = 0,
      sourceSink = T
      ),
    compoundFactory(
      compoundName = "N2",
      molarRatios = c(N=2),
      initialMols = 0
      ),
    compoundFactory(
      compoundName = "Ox",
      molarRatios = c(O=1),
      initialMols = 0,
      sourceSink = T
      )
  )

myGangstas =
  c(
    myGangstas,
    processFactory(
      myGangstas,
      name = "AutAssimCO2",
      energyTerm = -3.5E-02,
      fromCompoundNames = list(C = "CO2", O = "CO2", O = "CO2"),
      toCompoundNames = list(C = "Aut", O = "Aut", O = "Ox"),
      molarTerms = list(C = 1, O = BioStoichO, O = 2 - BioStoichO),
      organismName = "Aut"
    ),
    processFactory(
      myGangstas,
      name = "MetAssimCH4",
      energyTerm = -1.09E-03,
      fromCompoundNames = list(C = "CH4", O = "Ox"),
      toCompoundNames = list(C = "Met", O = "Met"),
      molarTerms = list(C = 1, O = BioStoichO),
      organismName = "Met"
    ),
    processFactory(
      myGangstas,
      name = "AutAssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3", O = "NO3"),
      toCompoundNames = list(N = "Aut", O = "Ox"),
      molarTerms = list(N = 1, O = 3),
      organismName = "Aut"
    ),
    processFactory(
      myGangstas,
      name = "MetAssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3", O = "NO3"),
      toCompoundNames = list(N = "Met", O = "Ox"),
      molarTerms = list(N = 1, O = 3),
      organismName = "Met"
    ),
    processFactory(
      myGangstas,
      name = "AutAssimNH4",
      energyTerm = -3.18E-05,
      fromCompoundNames = list(N = "NH4"),
      toCompoundNames = list(N = "Aut"),
      molarTerms = list(N = 1),
      organismName = "Aut"
    ),
    processFactory(
      myGangstas,
      name = "MetAssimNH4",
      energyTerm = -3.18E-05,
      fromCompoundNames = list(N = "NH4"),
      toCompoundNames = list(N = "Met"),
      molarTerms = list(N = 1),
      organismName = "Met"
    ),
    processFactory(
      myGangstas,
      name = "AutNitrif",
      energyTerm = 3.485E-4, #((1.83E-04 *3) + 1.48E-04)/2,
      fromCompoundNames = list(N = "NH4", O = "O2", O = "O2"),
      toCompoundNames = list(N = "NO3", O = "NO3", O = "Ox"),
      molarTerms = list(N = 1, O = 3, O = 2),
      organismName = "Aut"
    ),
    processFactory(
      myGangstas,
      name = "MetMethaneOxid",
      energyTerm = 8.18E-4, #4.09E-04 * 2,
      fromCompoundNames = list(C = "CH4", O = "O2", O = "O2"),
      toCompoundNames = list(C = "CO2", O = "CO2", O = "Ox"),
      molarTerms = list(C = 1, O = 2, O = 2),
      organismName = "Met"
    ),
    processFactory(
      myGangstas,
      name = "AutDecay",
      energyTerm = 0,
      fromCompoundNames = list(C = "Aut", N = "Aut", O = "Aut"),
      toCompoundNames = list(C = "DOM", N = "DOM", O = "DOM"),
      molarTerms = list(C = 1, N = BioStoichN, O = BioStoichO),
      organismName = c("Aut")
    ),
    processFactory(
      myGangstas,
      name = "MetDecay",
      energyTerm = 0,
      fromCompoundNames = list(C = "Met", N = "Met", O = "Met"),
      toCompoundNames = list(C = "DOM", N = "DOM", O = "DOM"),
      molarTerms = list(C = 1, N = BioStoichN, O = BioStoichO),
      organismName = c("Met")
    )
  )

myEquations <- makeEquations(
  gangstaObjects = myGangstas
)

writeGangstaModel(
  equations = myEquations
)


##

lp.gangsta = readGangsta.lp()

lpStatus = suppressWarnings(lpSolveAPI::solve.lpExtPtr(lp.gangsta))

compoundDifDF = compoundDifs(myGangstas, lp.gangsta, simple = F)
poolDifDF = poolDifs(myGangstas, lp.gangsta, simple = FALSE)


###  This is a silly comment that just says nothing.  MASTER FIX. Ann Marie

gangstaTest = function() {

  ### NEED TO ALLOW "." to access the molarRatio of the source for massTerm of transformation.
  ### Another comment to change.
  compoundParams = list(
    list(compoundName = "Het", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -1),
    list(compoundName = "Aut", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2),
    list(compoundName = "Met", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -3),
    list(compoundName = "DOM", molarRatios = c(C=1, N=6/106), initialMols = 0),
    list(compoundName = "CH4", molarRatios = c(C=1), initialMols = 0),
    list(compoundName = "NH4", molarRatios = c(N=1), initialMols = 0),
    list(compoundName = "NO3", molarRatios = c(N=1), initialMols = 0),
    list(compoundName = "NO2", molarRatios = c(N=1), initialMols = 0),
    list(compoundName = "N2O", molarRatios = c(N=1), initialMols = 0),
    list(compoundName = "O2" , molarRatios = c(O=1), initialMols = 0),
    list(compoundName = "SO4", molarRatios = c(S=1), initialMols = 0),
    list(compoundName = "CO2", molarRatios = c(C=1), initialMols = 0, sourceSink = T),
    list(compoundName = "N2" , molarRatios = c(N=1), initialMols = 0, sourceSink = T),
    list(compoundName = "HS" , molarRatios = c(S=1), initialMols = 0, sourceSink = T),
    list(compoundName = "Ox" , molarRatios = c(O=1), initialMols = 0, sourceSink = T)
  )

  compounds = unlist(lapply(compoundParams, do.call, what = "compoundFactory"), recursive = F)

  processParams = list(
    list(
      gangstaObjects = compounds,
      processName = "AssimDOM",
      energyTerm = -4.32E-04,
      fromCompoundNames = list(C = "DOM", N = "DOM"),
      toCompoundNames = list(C = ".", N = "."),
      massTerms = list(C = 1, N = 6/106),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "AssimCO2",
      energyTerm = -3.5E-02,
      fromCompoundNames = list(C = "CO2"),
      toCompoundNames = list(C = "."),
      massTerms = list(C = 1),
      organismNames = "Aut"
    ),
    list(
      gangstaObjects = compounds,
      processName = "AssimCH4",
      energyTerm = -1.09E-03,
      fromCompoundNames = list(C = "CH4"),
      toCompoundNames = list(C = "."),
      massTerms = list(C = 1),
      organismNames = "Met"
    ),
    list(
      gangstaObjects = compounds,
      processName = "AssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3"),
      toCompoundNames = list(N = "."),
      massTerms = list(N = 1),
      organismNames = c("Het", "Aut", "Met")
    ),
    list(
      gangstaObjects = compounds,
      processName = "AssimNO2",
      energyTerm = -1.25E-04,
      fromCompoundNames = list(N = "NO2"),
      toCompoundNames = list(N = "."),
      massTerms = list(N = 1),
      organismNames = c("Het", "Aut", "Met")
    ),
    list(
      gangstaObjects = compounds,
      processName = "AssimNH4",
      energyTerm = -3.18E-05,
      fromCompoundNames = list(N = "NH4"),
      toCompoundNames = list(N = "."),
      massTerms = list(N = 1),
      organismNames = c("Het", "Aut", "Met"),
      limitToInitMass = c(F, T, T)
    ),
    list(
      gangstaObjects = compounds,
      processName = "AerobicOfDOM",
      energyTerm = 4.37-04,
      fromCompoundNames = list(C = "DOM", N = "DOM", O = "O2"),
      toCompoundNames = list(C = "CO2", N = "NH4", O = "Ox"),
      massTerms = list(C = 1, N = 6/106, O = 2),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "AerobicOfHet",
      energyTerm = 4.37-04,
      fromCompoundNames = list(C = ".", N = ".", O = "O2"),
      toCompoundNames = list(C = "CO2", N = "NH4", O = "Ox"),
      massTerms = list(C = 1, N = 16/106, O = 2),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "DenitOfHetStep1",
      energyTerm = 2.88E-04,
      fromCompoundNames = list(C = ".", N = ".", N = "NO3"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "NO2"),
      massTerms = list(C = 1, N = 16/106, N = 2),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "DenitOfDOMStep1",
      energyTerm = 2.88E-04,
      fromCompoundNames = list(C = "DOM", N = "DOM", N = "NO3"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "NO2"),
      massTerms = list(C = 1, N = 16/106, N = 2),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "DenitOfHetStep2",
      energyTerm = 4.15E-04,
      fromCompoundNames = list(C = ".", N = ".", N = "NO2"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2O"),
      massTerms = list(C = 1, N = 16/106, N = 2),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "DenitOfDOMStep2",
      energyTerm = 4.15E-04,
      fromCompoundNames = list(C = "DOM", N = "DOM", N = "NO2"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2O"),
      massTerms = list(C = 1, N = 16/106, N = 2),
      organismNames = "Het"
    ),

    list(
      gangstaObjects = compounds,
      processName = "DenitOfHetStep3",
      energyTerm = 6.45E-04,
      fromCompoundNames = list(C = ".", N = ".", N = "N2O"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2"),
      massTerms = list(C = 1, N = 16/106, N = 4),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "DenitOfDOMStep3",
      energyTerm = 6.45E-04,
      fromCompoundNames = list(C = "DOM", N = "DOM", N = "N2O"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2"),
      massTerms = list(C = 1, N = 16/106, N = 4),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "SulfateRedOfHet",
      energyTerm = 3.8E-05,
      fromCompoundNames = list(C = ".", N = ".", S = "SO4"),
      toCompoundNames = list(C = "CO2", N = "NH4", S = "HS"),
      massTerms = list(C = 1, N = 16/106, S = 0.5),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "SulfateRedOfDOM",
      energyTerm = 3.8E-05,
      fromCompoundNames = list(C = "DOM", N = "DOM", S = "SO4"),
      toCompoundNames = list(C = "CO2", N = "NH4", S = "HS"),
      massTerms = list(C = 1, N = 16/106, S = 0.5),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "MethanogenesisOfHet",
      energyTerm = 2.8E-05,
      fromCompoundNames = list(C = ".", C = ".", N = "."),
      toCompoundNames = list(C = "CO2", C = "CH4", N = "NH4"),
      massTerms = list(C = 0.5, C = 0.5, N = 16/106),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "MethanogenesisOfDOM",
      energyTerm = 2.8E-05,
      fromCompoundNames = list(C = "DOM", C = "DOM", N = "DOM"),
      toCompoundNames = list(C = "CO2", C = "CH4", N = "NH4"),
      massTerms = list(C = 0.5, C = 0.5, N = 16/106),
      organismNames = "Het"
    ),
    list(
      gangstaObjects = compounds,
      processName = "NitrifStep1",
      energyTerm = 1.83E-04,
      fromCompoundNames = list(N = "NH4", O = "O2"),
      toCompoundNames = list(N = "NO2", O = "Ox"),
      massTerms = list(N = 2/3, O = 2),
      organismNames = "Aut"
    ),
    list(
      gangstaObjects = compounds,
      processName = "NitrifStep2",
      energyTerm = 1.48E-04,
      fromCompoundNames = list(N = "NO2", O = "O2"),
      toCompoundNames = list(N = "NO3", O = "Ox"),
      massTerms = list(N = 2, O = 2),
      organismNames = "Aut"
    ),
    list(
      gangstaObjects = compounds,
      processName = "MethaneOxid",
      energyTerm = 4.09E-04,
      fromCompoundNames = list(C = "CH4", O = "O2"),
      toCompoundNames = list(C = "CO2", O = "Ox"),
      massTerms = list(C = 0.5, O = 2),
      organismNames = "Met"
    )
  )

  processes = unlist(lapply(processParams, do.call, what = processFactory), recursive = F)
  return(c(compounds, processes))
}



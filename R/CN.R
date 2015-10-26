CNModel = function() {
  compoundParams = list(
    list(compoundName = "Het", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "Aut", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "Met", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "DOM", molarRatios = c(C=1, N=16/106), initialMols = 0),
    list(compoundName = "CH4", molarRatios = c(C=1), initialMols = 0),
    list(compoundName = "NH4", molarRatios = c(N=1), initialMols = 0),
    list(compoundName = "NO3", molarRatios = c(N=1), initialMols = 0),
#    list(compoundName = "O2" , molarRatios = c(O=1), initialMols = 0, sourceSink = T),
    list(compoundName = "CO2", molarRatios = c(C=1), initialMols = 0, sourceSink = T),
    list(compoundName = "N2" , molarRatios = c(N=1), initialMols = 0, sourceSink = T)
#    list(compoundName = "Ox" , molarRatios = c(O=1), initialMols = 0, sourceSink = T)
  )

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
      energyTerm = 4.37E-04,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM")), #, O = "O2"),
      toCompoundNames = list(C = "CO2", N = "NH4"), #, O = "Ox"),
      molarTerms = list(C = 1, N = NA), #, O = 2),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "Denit",
      energyTerm = (2 * 2.88E-04) + (2 * 4.15E-04) + 6.45E-04,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "NO3"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2"),
      molarTerms = list(C = 5, N = NA, N = 4),
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
      energyTerm = (1.83E-04 *3) + 1.48E-04,
      fromCompoundNames = list(N = "NH4"), #, O = "O2"),
      toCompoundNames = list(N = "NO3"), #, O = "Ox"),
      molarTerms = list(N = 2), #, O = 8),
      organismName = "Aut"
    ),
    list(
      name = "MethaneOxid",
      energyTerm = 4.09E-04,
      fromCompoundNames = list(C = "CH4"), #, O = "O2"),
      toCompoundNames = list(C = "CO2"), #, O = "Ox"),
      molarTerms = list(C = 0.5), #, O = 2),
      organismName = "Met"
    ),
    list(
      name = "Decay",
      energyTerm = 0,
      fromCompoundNames = list(C = ".", N = "."),
      toCompoundNames = list(C = "DOM", N = "DOM"),
      molarTerms = list(C = 1, N = NA),
      organismName = c("Aut", "Met")
    )
  )

  tag = "CN"
  doAll(tag, compoundParams, processParams, O = F, S = F)
}


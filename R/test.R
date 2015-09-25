testList = function() {
  return(
    list(
      name = "Aerobic",
      energyTerm = 4.37-04,
      fromCompoundNames = list(C = c("DOM", "."), N = c("DOM", "."), O = "O2"),
      toCompoundNames = list(C = "CO2", N = "NH4", O = "Ox"),
      organismName = c("Het", "Aut", "Met"),
      molarTerms = list(C = 1, N = NA, O = 2),
      limitToInitMols = c(F, T, T),
      processSuffix = c("ofDOM", "ofHet")
    )
  )
}

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
      name = "AssimNO2",
      energyTerm = -1.25E-04,
      fromCompoundNames = list(N = "NO2"),
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
      energyTerm = 4.37-04,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), O = "O2"),
      toCompoundNames = list(C = "CO2", N = "NH4", O = "Ox"),
      molarTerms = list(C = 1, N = NA, O = 2),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "DenitStep1",
      energyTerm = 2.88E-04,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "NO3"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "NO2"),
      molarTerms = list(C = 1, N = NA, N = 2),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "DenitStep2",
      energyTerm = 4.15E-04,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "NO2"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2O"),
      molarTerms = list(C = 1, N = NA, N = 2),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "DenitStep3",
      energyTerm = 6.45E-04,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "N2O"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2"),
      molarTerms = list(C = 1, N = NA, N = 4),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "SulfateRed",
      energyTerm = 3.8E-05,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), S = "SO4"),
      toCompoundNames = list(C = "CO2", N = "NH4", S = "HS"),
      molarTerms = list(C = 1, N = NA, S = 0.5),
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
      name = "NitrifStep1",
      energyTerm = 1.83E-04,
      fromCompoundNames = list(N = "NH4", O = "O2"),
      toCompoundNames = list(N = "NO2", O = "Ox"),
      molarTerms = list(N = 2/3, O = 2),
      organismName = "Aut"
    ),
    list(
      name = "NitrifStep2",
      energyTerm = 1.48E-04,
      fromCompoundNames = list(N = "NO2", O = "O2"),
      toCompoundNames = list(N = "NO3", O = "Ox"),
      molarTerms = list(N = 2, O = 2),
      organismName = "Aut"
    ),
    list(
      name = "MethaneOxid",
      energyTerm = 4.09E-04,
      fromCompoundNames = list(C = "CH4", O = "O2"),
      toCompoundNames = list(C = "CO2", O = "Ox"),
      molarTerms = list(C = 0.5, O = 2),
      organismName = "Met"
    )
  )

#   processParams = list(
#     list(
#       name = "Methanogenesis",
#       energyTerm = 2.8E-05,
#       fromCompoundNames = list(C = c(".", "DOM"), C = c(".", "."), N = c(".", "DOM"), N = c("NH4", ".")),
#       toCompoundNames = list(C = "CO2", C = "CH4", N = "NH4", N = c("NO3", "NH4")),
#       molarTerms = list(C = 0.5, C = 0.5, N = c(NA, NA), N = c(1, NA)),
#       organismName = "Het",
#       processSuffix = c("ofHet", "ofDOM")
#     )
#   )

  expandedSpecs = unlist(lapply(processParams, expandMultiprocessSpec), recursive = F)
  processes = unlist(lapply(expandedSpecs, function(x) do.call(processFactory, args = c(list(compounds), x))), recursive = F)
  names(processes) = getGangstaAttribute(processes, gangstaAttributeName("name"))
  return(c(compounds, processes))
}



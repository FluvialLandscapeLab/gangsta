CNOS_Any = function(activeElements, sourceSinks = c("CO2", "N2", "HS", "Ox")) {

  tag = paste(activeElements, collapse = "")

  compoundParams = list(
    list(compoundName = "Het", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "Aut", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "Met", molarRatios = c(C=1, N=16/106), initialMols = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "DOM", molarRatios = c(C=1, N=16/106), initialMols = 0),
    list(compoundName = "CH4", molarRatios = c(C=1),           initialMols = 0),
    list(compoundName = "NH4", molarRatios = c(N=1),           initialMols = 0),
    list(compoundName = "NO3", molarRatios = c(N=1, O=3),      initialMols = 0),
    list(compoundName = "O2" , molarRatios = c(O=1),           initialMols = 0),
    list(compoundName = "SO4", molarRatios = c(S=1, O=4),      initialMols = 0),
    list(compoundName = "CO2", molarRatios = c(C=1, O=2),      initialMols = 0),
    list(compoundName = "N2" , molarRatios = c(N=1),           initialMols = 0),
    list(compoundName = "HS" , molarRatios = c(S=1),           initialMols = 0),
    list(compoundName = "Ox" , molarRatios = c(O=1),           initialMols = 0)
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
      fromCompoundNames = list(C = "CO2", O = "CO2"),
      toCompoundNames = list(C = ".", O = "Ox"),
      molarTerms = list(C = 1, O = 2),
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
      fromCompoundNames = list(N = "NO3", O = "NO3"),
      toCompoundNames = list(N = ".", O = "Ox"),
      molarTerms = list(N = 1, O = 3),
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
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), O = "O2"),
      toCompoundNames = list(C = "CO2", N = "NH4", O = "CO2"),
      molarTerms = list(C = 1, N = NA, O = 2),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "Denit",
      energyTerm = ((2 * 2.88E-04) + (2 * 4.15E-04) + 6.45E-04)/5,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), N = "NO3", O = "NO3", O = "NO3"),
      toCompoundNames = list(C = "CO2", N = "NH4", N = "N2", O = "CO2", O = "Ox"),
      molarTerms = list(C = 1, N = NA, N = 4/5, O = 2, O = 2/5),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "SulfateRed",
      energyTerm = 3.8E-05,
      fromCompoundNames = list(C = c(".", "DOM"), N = c(".", "DOM"), S = "SO4", O = "SO4"),
      toCompoundNames = list(C = "CO2", N = "NH4", S = "HS", O = "CO2"),
      molarTerms = list(C = 1, N = NA, S = 0.5, O = 2),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "Methanogenesis",
      energyTerm = 2.8E-05,
      fromCompoundNames = list(C = c(".", "DOM"), C = c(".", "DOM"), N = c(".", "DOM"), O = "Ox"),
      toCompoundNames = list(C = "CO2", C = "CH4", N = "NH4", O = "CO2"),
      molarTerms = list(C = 0.5, C = 0.5, N = NA, O = 1),
      organismName = "Het",
      processSuffix = c("ofHet", "ofDOM")
    ),
    list(
      name = "Nitrif",
      energyTerm = ((1.83E-04 *3) + 1.48E-04)/2,
      fromCompoundNames = list(N = "NH4", O = "O2", O = "O2"),
      toCompoundNames = list(N = "NO3", O = "NO3", O = "Ox"),
      molarTerms = list(N = 1, O = 3, O = 1),
      organismName = "Aut"
    ),
    list(
      name = "MethaneOxid",
      energyTerm = 4.09E-04 * 2,
      fromCompoundNames = list(C = "CH4", O = "O2", O = "O2"),
      toCompoundNames = list(C = "CO2", O = "CO2", O = "Ox"),
      molarTerms = list(C = 1, O = 2, O = 2),
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

  keep = sapply(compoundParams, function(x) names(x[["molarRatios"]])[1] %in% activeElements)
  compoundParams = compoundParams[keep]
  compoundParams = lapply(
    compoundParams,
    function(x) {
      x[["molarRatios"]] = x[["molarRatios"]][names(x[["molarRatios"]]) %in% activeElements]
      x[["molarRatios"]] = x[["molarRatios"]] / x[["molarRatios"]][1] # set molar ratio of first element to 1.0
      if (x[["compoundName"]] %in% sourceSinks) x = c(x, list(sourceSink = T))
      return(x)
    }
  )

  compoundNames = sapply(compoundParams, "[[", "compoundName")

  processParams = unlist(lapply(processParams, expandMultiprocessSpec), recursive = F)

  keep = sapply(processParams, function(x) any(names(x[["fromCompoundNames"]]) %in% activeElements))
  processParams = processParams[keep]
  kill = sapply(
    processParams,
    function(x) {
      froms = x[["fromCompoundNames"]]
      froms[froms == "."] = x[["organismName"]]
      noGood = any(!(froms %in% compoundNames) & !(froms %in% sourceSinks))
      return(noGood)
    }
  )
  processParams = processParams[!kill]

  processParams = lapply(
    processParams,
    function(x) {
      for(i in c("fromCompoundNames", "toCompoundNames", "molarTerms")) {
        x[[i]] = x[[i]][names(x[[i]]) %in% activeElements]
      }
      return(x)
    }
  )

  doAll(tag, compoundParams, processParams, compoundNames)
}


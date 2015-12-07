CNOSH_Any = function(activeElements, sourceSinks = character(0)) {

  sourceSinks = c("CO2", "Ox", "Hx", sourceSinks)

  tag = paste(activeElements, collapse = "")

  BioStoichN = 16/106
  BioStoichO = 110/106
  BioStoichS = 1.7/106
  BioStoichH = 263/106
  DOMStoichN = 16/106
  DOMStoichO = 110/106
  DOMStoichS = 1.7/106
  DOMStoichH = 263/106

  compoundParams = list(
    list(compoundName = "Het" , molarRatios = c(C=1, N=BioStoichN, O=BioStoichO, S=BioStoichS, H = BioStoichH), initialMols = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "Aut" , molarRatios = c(C=1, N=BioStoichN, O=BioStoichO, S=BioStoichS, H = BioStoichH), initialMols = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "Met" , molarRatios = c(C=1, N=BioStoichN, O=BioStoichO, S=BioStoichS, H = BioStoichH), initialMols = 0, respirationRate = -2.83E-6 * 24),
    #    list(compoundName = "DOMX", molarRatios = c(C=1, N=BioStoichN, O=BioStoichO), initialMols = 0),
    list(compoundName = "DOM" , molarRatios = c(C=1, N=DOMStoichN, O=DOMStoichO, S=DOMStoichS, H = DOMStoichH), initialMols = 0),
    list(compoundName = "CH4" , molarRatios = c(C=1, H=4),      initialMols = 0),
    list(compoundName = "NH4" , molarRatios = c(N=1, H=4),      initialMols = 0),
    list(compoundName = "NO3" , molarRatios = c(N=1, O=3),      initialMols = 0),
    list(compoundName = "O2"  , molarRatios = c(O=2),           initialMols = 0),
    list(compoundName = "SO4" , molarRatios = c(S=1, O=4),      initialMols = 0),
    list(compoundName = "CO2" , molarRatios = c(C=1, O=2),      initialMols = 0),
    list(compoundName = "N2"  , molarRatios = c(N=2),           initialMols = 0),
    list(compoundName = "HS"  , molarRatios = c(S=1, H=1),      initialMols = 0),
    list(compoundName = "Ox"  , molarRatios = c(O=1),           initialMols = 0),
    list(compoundName = "Hx"  , molarRatios = c(H=1),           initialMols = 0)
  )

  processParams = list(
    list(
      name = "AssimDOM",
      energyTerm = -4.32E-04,
      fromCompoundNames = list(
        C = "DOM",
        N = "DOM",
        O = "DOM",
        S = "DOM",
        H = "DOM"),
      toCompoundNames = list(
        C = ".",
        N = ".",
        O = ".",
        S = ".",
        H = "."),
      molarTerms = list(
        C = 1,
        N = DOMStoichN,
        O = DOMStoichO,
        S = DOMStoichS,
        H = DOMStoichH),
      organismName = "Het"
    ),
    list(
      name = "AssimCO2",
      energyTerm = -3.5E-02,
      fromCompoundNames = list(C = "CO2", O = "CO2", O = "CO2", H = "Hx"),
      toCompoundNames = list(C = ".", O = ".", O = "Ox", H = "."),
      molarTerms = list(C = 1, O = BioStoichO, O = 2 - BioStoichO, H = BioStoichH),
      organismName = "Aut"
    ),
    list(
      name = "AssimCH4",
      energyTerm = -1.09E-03,
      fromCompoundNames = list(C = "CH4", H = "CH4", H = "CH4", O = "Ox"),
      toCompoundNames = list(C = ".", H = ".", H = "Hx", O = "."),
      molarTerms = list(C = 1, H = BioStoichH, H = 4 - BioStoichH, O = BioStoichO),
      organismName = "Met"
    ),
    list(
      name = "AssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3", O = "NO3"),
      # O must go to Ox, or else in an Oxygen only model, microbes can grow by ammimilating NO3
      toCompoundNames = list(N = ".", O = "Ox"),
      molarTerms = list(N = 1, O = 3),
      organismName = c("Het", "Aut", "Met")
    ),
    list(
      name = "AssimNH4",
      energyTerm = -3.18E-05,
      fromCompoundNames = list(N = "NH4", H = "NH4"),
      # H must go to Hx, or else in an Hydrogen only model, microbes can grow by assimilating NH4
      toCompoundNames = list(N = ".", H = "Hx"),
      molarTerms = list(N = 1, H = 4),
      organismName = c("Het", "Aut", "Met"),
      limitToInitMols = c(F, T, T)
    ),
    list(
      name = "AssimSO4",
      energyTerm = -1.5E-04,  ## Value is made up -- need a real value
      fromCompoundNames = list(S = "SO4", O = "SO4"),
      # O must go to Ox, or else in an Oxygen only model, microbes can grow by ammimilating SO4
      toCompoundNames = list(S = ".", O = "Ox"),
      molarTerms = list(S = 1, O = 4),
      organismName = c("Het", "Aut", "Met")
    ),
    list(
      name = "AssimHS",
      energyTerm = -3E-05,  ## Value is made up -- need a real value
      fromCompoundNames = list(S = "HS", H = "HS"),
      # H must go to Hx, or else in an Hydrogen only model, microbes can grow by ammimilating HS
      toCompoundNames = list(S = ".", H = "Hx"),
      molarTerms = list(S = 1, H = 1),
      organismName = c("Het", "Aut", "Met")
    ),
    list(
      name = "Aerobic",
      energyTerm = 4.37E-04,
      fromCompoundNames = list(
        C = c(".", "DOM"),
        N = c(".", "DOM"),
        O = c(".", "DOM"),
        S = c(".", "DOM"),
        H = c(".", "DOM"),

        H = "Hx",
        H = "Hx",
        O = "O2"),
      toCompoundNames = list(
        C = "CO2",
        N = "NH4",
        O = "Ox",
        S = "HS",
        H = "Hx",

        H = "HS",
        H = "NH4",
        O = "CO2"),
      molarTerms = list(
        C = 1,
        N = c(BioStoichN, DOMStoichN),
        O = c(BioStoichO, DOMStoichO),
        S = c(BioStoichS, DOMStoichS),
        H = c(BioStoichH, DOMStoichH),

        H = c(BioStoichS, DOMStoichS),
        H = c(BioStoichN * 4, DOMStoichN * 4),
        O = 2),
      organismName = "Het",
      processSuffix = c("ofBiomass", "ofDOM")
    ),
    list(
      name = "Denit",
      energyTerm = 4.102E-4, # ((2 * 2.88E-04) + (2 * 4.15E-04) + 6.45E-04)/5,
      fromCompoundNames = list(
        C = c(".", "DOM"),
        N = c(".", "DOM"),
        O = c(".", "DOM"),
        S = c(".", "DOM"),
        H = c(".", "DOM"),

        N = "NO3",
        O = "NO3",
        O = "NO3",
        H = "Hx",
        H = "Hx"),
      toCompoundNames = list(
        C = "CO2",
        N = "NH4",
        O = "Ox",
        S = "HS",
        H = "Hx",

        N = "N2",
        O = "CO2",
        O = "Ox",
        H = "NH4",
        H = "HS"),
      molarTerms = list(
        C = 1,
        N = c(BioStoichN, DOMStoichN),
        O = c(BioStoichO, DOMStoichO),
        S = c(BioStoichS, DOMStoichS),
        H = c(BioStoichH, DOMStoichH),

        N = 4/5,
        O = 2,
        O = 2/5,
        H = c(BioStoichN * 4, DOMStoichN * 4),
        H = c(BioStoichS, DOMStoichS)
      ),
      organismName = "Het",
      processSuffix = c("ofBiomass", "ofDOM")
    ),
    list(
      name = "SulfateRed",
      energyTerm = 3.8E-05,
      fromCompoundNames = list(
        C = c(".", "DOM"),
        N = c(".", "DOM"),
        O = c(".", "DOM"),
        S = c(".", "DOM"),
        H = c(".", "DOM"),

        H = "Hx",

        S = "SO4",
        O = "SO4",
        H = "Hx"),
      toCompoundNames = list(
        C = "CO2",
        N = "NH4",
        O = "Ox",
        S = "HS",
        H = "Hx",

        H = "NH4",

        S = "HS",
        O = "CO2",
        H = "HS"),
      molarTerms = list(
        C = 1,
        N = c(BioStoichN, DOMStoichN),
        O = c(BioStoichO, DOMStoichO),
        S = c(BioStoichS, DOMStoichS),
        H = c(BioStoichH, DOMStoichH),

        H = c(BioStoichN * 4, DOMStoichN * 4),

        S = 0.5,
        O = 2,
        H = 0.5 + c(BioStoichS, DOMStoichS)),
      organismName = "Het",
      processSuffix = c("ofBiomass", "ofDOM")
    ),
    list(
      name = "Methanogenesis",
      energyTerm = 2.8E-05,
      fromCompoundNames = list(
        C = c(".", "DOM"),
        C = c(".", "DOM"),
        N = c(".", "DOM"),
        O = c(".", "DOM"),
        S = c(".", "DOM"),
        H = c(".", "DOM"),

        O = "Ox",
        H = "Hx",
        H = "Hx",
        H = "Hx"),
      toCompoundNames = list(
        C = "CO2",
        C = "CH4",
        N = "NH4",
        O = "Ox",
        S = "HS",
        H = "Hx",

        O = "CO2",
        H = "CH4",
        H = "NH4",
        H = "HS"),
      molarTerms = list(
        C = 0.5,
        C = 0.5,
        N = c(BioStoichN, DOMStoichN),
        O = c(BioStoichO, DOMStoichO),
        S = c(BioStoichS, DOMStoichS),
        H = c(BioStoichH, DOMStoichH),

        O = 1,
        H = 2,
        H = c(BioStoichN * 4, DOMStoichN * 4),
        H = c(BioStoichS, BioStoichS)),
      organismName = "Het",
      processSuffix = c("ofBiomass", "ofDOM")
    ),
    list(
      name = "Nitrif",
      energyTerm = 3.485E-4, #((1.83E-04 *3) + 1.48E-04)/2,
      fromCompoundNames = list(N = "NH4", H = "NH4", O = "O2", O = "O2"),
      toCompoundNames = list(N = "NO3", H = "Hx", O = "NO3", O = "Ox"),
      molarTerms = list(N = 1, H = 4, O = 3, O = 1),
      organismName = "Aut"
    ),
    list(
      name = "MethaneOxid",
      energyTerm = 8.18E-4, #4.09E-04 * 2,
      fromCompoundNames = list(C = "CH4", H = "CH4", O = "O2", O = "O2"),
      toCompoundNames = list(C = "CO2", H = "Hx", O = "CO2", O = "Ox"),
      molarTerms = list(C = 1, H = 4, O = 2, O = 2),
      organismName = "Met"
    ),
    list(
      name = "Decay",
      energyTerm = 0,
      fromCompoundNames = list(C = ".", N = ".", O = ".", S = ".", H = "."),
      toCompoundNames = list(C = "DOM", N = "DOM", O = "DOM", S = "DOM", H = "DOM"),
      molarTerms = list(C = 1, N = BioStoichN, O = BioStoichO, S = BioStoichS, H = BioStoichH),
      organismName = c("Het", "Aut", "Met")
    )
  )

  keep = sapply(compoundParams, function(x) any(names(x[["molarRatios"]]) %in% activeElements))
  compoundParams = compoundParams[keep]
  compoundParams = lapply(
    compoundParams,
    function(x) {
      x[["molarRatios"]] = x[["molarRatios"]][names(x[["molarRatios"]]) %in% activeElements]
      # set molar ratio of first remaining element (reference element) to 1.0
      # x[["molarRatios"]] = x[["molarRatios"]] / x[["molarRatios"]][1]
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

  doAll(tag, compoundParams, processParams, compoundNames, sourceSinks)
}


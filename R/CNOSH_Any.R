CNOSH_Any = function(activeElements, InfiniteCompounds = c("Ox", "Hx")) {

  activeElementTag = paste(activeElements, collapse = "")
  InfiniteCompoundTag = paste(InfiniteCompounds, collapse = ".")

  modelNameTag = paste0(activeElementTag, "_", InfiniteCompoundTag, collapse = "")

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
    list(compoundName = "Het" , molarRatios = c(C=1, N=BioStoichN, O=BioStoichO, S=BioStoichS, H = BioStoichH), initialMolecules = 0, respirationRate = -2.83E-6 * 24), #respirationRate units kJ umol-1 Day-1
    list(compoundName = "Aut" , molarRatios = c(C=1, N=BioStoichN, O=BioStoichO, S=BioStoichS, H = BioStoichH), initialMolecules = 0, respirationRate = -2.83E-6 * 24),
    list(compoundName = "Met" , molarRatios = c(C=1, N=BioStoichN, O=BioStoichO, S=BioStoichS, H = BioStoichH), initialMolecules = 0, respirationRate = -2.83E-6 * 24),
    #    list(compoundName = "DOMX", molarRatios = c(C=1, N=BioStoichN, O=BioStoichO), initialMolecules = 0),
    list(compoundName = "DOM" , molarRatios = c(C=1, N=DOMStoichN, O=DOMStoichO, S=DOMStoichS, H = DOMStoichH), initialMolecules = 0),
    list(compoundName = "CH4" , molarRatios = c(C=1, H=4),      initialMolecules = 0),
    list(compoundName = "NH4" , molarRatios = c(N=1, H=4),      initialMolecules = 0),
    list(compoundName = "NO3" , molarRatios = c(N=1, O=3),      initialMolecules = 0),
    list(compoundName = "O2"  , molarRatios = c(O=2),           initialMolecules = 0),
    list(compoundName = "SO4" , molarRatios = c(S=1, O=4),      initialMolecules = 0),
    list(compoundName = "CO2" , molarRatios = c(C=1, O=2),      initialMolecules = 0),
    list(compoundName = "N2"  , molarRatios = c(N=2),           initialMolecules = 0),
    list(compoundName = "HS"  , molarRatios = c(S=1, H=1),      initialMolecules = 0),
    list(compoundName = "Ox"  , molarRatios = c(O=1),           initialMolecules = 0),
    list(compoundName = "Hx"  , molarRatios = c(H=1),           initialMolecules = 0)
  )

  processParams = list(
    list(
      name = "AssimO",
      energyTerm = 0,
      fromCompoundNames = list(O = "Ox"),
      toCompoundNames = list(O = "."),
      molarTerms = list(O = 1),
      organismName = c("Het", "Aut", "Met")
    ),
    list(
      name = "AssimH",
      energyTerm = 0,
      fromCompoundNames = list(H = "Hx"),
      toCompoundNames = list(H = "."),
      molarTerms = list(H = 1),
      organismName = c("Het", "Aut", "Met")
    ),
    list(
      name = "AssimDOM",
      energyTerm = -4.32E-04, # units are kJ *(umols of compound)-1
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
      fromCompoundNames = list(C = "CO2", O = "CO2"),
      toCompoundNames = list(C = ".", O = c(".", "Ox")),
      molarTerms = list(C = 1, O = 2),
      organismName = "Aut"
    ),
    list(
      name = "AssimCH4",
      energyTerm = -1.09E-03,
      fromCompoundNames = list(C = "CH4", H = "CH4"),
      toCompoundNames = list(C = ".", H = c(".", "Hx")),
      molarTerms = list(C = 1, H = 4),
      organismName = "Met"
    ),
    list(
      name = "AssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3", O = "NO3"),
      toCompoundNames = list(N = ".", O = c(".", "Ox")),
      molarTerms = list(N = 1, O = 3),
      organismName = c("Het", "Aut", "Met")
    ),
    list(
      name = "AssimNH4",
      energyTerm = -3.18E-05,
      fromCompoundNames = list(N = "NH4", H = "NH4"),
      toCompoundNames = list(N = ".", H = c(".", "Hx")),
      molarTerms = list(N = 1, H = 4),
      organismName = c("Het", "Aut", "Met"),
      limitToInitMolecules = c(F, T, T)
    ),
    list(
      name = "AssimSO4",
      energyTerm = -9.28E-05,  ## Based on PAPS pathway in Shen and Buick (2004)
      fromCompoundNames = list(S = "SO4", O = "SO4"),
      toCompoundNames = list(S = ".", O = c(".", "Ox")),
      molarTerms = list(S = 1, O = 4),
      organismName = c("Het", "Aut", "Met")
    ),
    list(
      name = "AssimHS",
      energyTerm = -1E-7,  ## Based on Shen and Buick (2004) Fig 1, this could be zero
      fromCompoundNames = list(S = "HS", H = "HS"),
      toCompoundNames = list(S = ".", H = c(".", "Hx")),
      molarTerms = list(S = 1, H = 1),
      organismName = c("Het", "Aut", "Met"),
      limitToInitMolecules = c(F, T, T)
    ),
    list(
      name = "Aerobic",
      energyTerm = 4.37E-04,
      fromCompoundNames = list(
        C = list(".", "DOM"),
        N = list(".", "DOM"),
        O = list(".", "DOM"),
        S = list(".", "DOM"),
        H = list(".", "DOM"),

        H = "Hx",
        H = "Hx",
        O = "O2"
        ),
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
        C = list(".", "DOM"),
        N = list(".", "DOM"),
        O = list(".", "DOM"),
        S = list(".", "DOM"),
        H = list(".", "DOM"),

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
        C = list(".", "DOM"),
        N = list(".", "DOM"),
        O = list(".", "DOM"),
        S = list(".", "DOM"),
        H = list(".", "DOM"),

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
        H = 0.5 + c(BioStoichS, DOMStoichS)), # Hx goes to HS when we make half a mol of HS PLUS we need to account for the HS-H produced from the breakdown of OM
      organismName = "Het",
      processSuffix = c("ofBiomass", "ofDOM")
    ),
    list(
      name = "SulfideOxidation",
      energyTerm = 3.983E-04,
      fromCompoundNames = list(
        H = "HS",
        S = "HS",
        O = "O2"
      ),
      toCompoundNames = list(
        H = "Hx",
        S = "SO4",
        O = "SO4"
      ),
      molarTerms = list(
        H = 0.5,
        S = 0.5,
        O = 2
      ),
      organismName = "Aut"
    ),
    list(
      name = "Methanogenesis",
      energyTerm = 2.8E-05,
      fromCompoundNames = list(
        C = list(".", "DOM"),
        C = list(".", "DOM"),
        N = list(".", "DOM"),
        O = list(".", "DOM"),
        S = list(".", "DOM"),
        H = list(".", "DOM"),

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

  names(compoundParams) = sapply(compoundParams, "[[", "compoundName")

  # tag for inclusion any compoundParams that have at least one active element in the molarRatios vector
  keep = sapply(compoundParams, function(x) any(names(x[["molarRatios"]]) %in% activeElements))
  compoundParams = compoundParams[keep]

  #
  compoundParams = lapply(
    compoundParams,
    function(x) {
      # keep the molar ratios for active elements ONLY
      x[["molarRatios"]] = x[["molarRatios"]][names(x[["molarRatios"]]) %in% activeElements]
      # for any compound in the InfiniteCompounds list, add the InfiniteCompound = T attribute to the compoundParams list
      if (x[["compoundName"]] %in% InfiniteCompounds) x = c(x, list(InfiniteCompound = T))
      return(x)
    }
  )

  # make a vector of remaining compound names
  compoundNames = names(compoundParams)

  # expand multi-process params into process params.  In other words, each
  # process is replicated for each organismName vector and the organismName is
  # appended to the process name.  Also, replaces "." in from/to compound vector
  # with name of organism.
  processParams = unlist(lapply(processParams, expandMultiprocessSpec), recursive = F)

  # Tag for inclusion all processes where from compounds involve active
  # elements.
  keep = sapply(processParams, function(x) any(names(x[["fromCompoundNames"]]) %in% activeElements))
  processParams = processParams[keep]

  # Tag for exclusion any process where the from compound has been removed from
  # the compound list AND is not a source/sink.
  kill = sapply(
    processParams,
    function(x) {
      froms = x[["fromCompoundNames"]]
      froms[froms == "."] = x[["organismName"]]
      noGood = any(!(froms %in% compoundNames) & !(froms %in% InfiniteCompounds))
      return(noGood)
    }
  )
  processParams = processParams[!kill]

  # this local function is used in the lapply, below.  It changes any indexes in
  # the transfer groups to 0 if they are associated with transfers of elements
  # that are not active and subtracts 1 from all larger indexes.  keep is a
  # logical of whether the element of a transfer is active.
  keepTransferOptions = function (transferOptions, keep) {
    #get indexs of transfers to kill
    kill = which(!keep)
    #need to kill indexes from largest to smallest.  killing smaller index first
    #changes the value of any larger indexes that need to be removed.
    if(length(kill) > 0) kill = kill[length(kill):1]
    for(k in kill) {
      transferOptions =
        lapply (
          transferOptions,
          function(opt) {
            opt[opt == k] = 0
            opt[opt > k] = opt[opt > k] - 1
            return(opt)
          }
        )
    }
    transferOptions = lapply(transferOptions, function(opt) opt = opt[opt > 0])
    transferOptions = transferOptions[sapply(transferOptions, length) > 0]
    return(transferOptions)
  }

  # Remove any pools that contain elements not in the activeElements list.
  processParams =
    lapply(
      processParams,
      function(pp) {
        keep = names(pp$fromCompoundNames) %in% activeElements
        pp$fromCompoundNames = pp$fromCompoundNames[keep]
        pp$toCompoundNames = pp$toCompoundNames[keep]
        pp$molarTerms = pp$molarTerms[keep]
        pp$transferOptions = keepTransferOptions(pp$transferOptions, keep)
        return(pp)
      }
    )

  doAll(tag, modelNameTag, compoundParams, processParams, compoundNames, InfiniteCompounds)
}


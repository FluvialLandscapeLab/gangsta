organism = list(molarRatios = c(C=1, N=16/106), respirationRate = -99)

gangstaTest = function() {
  compoundParams = list(
    c(list(compoundName = "Het"), organism),
    c(list(compoundName = "Aut"), organism),
    c(list(compoundName = "Met"), organism),
    list(compoundName = "DOM", molarRatios = c(C=1, N=16/106)),
    list(compoundName = "CH4", molarRatios = c(C=1)),
    list(compoundName = "NH4", molarRatios = c(N=1)),
    list(compoundName = "NO3", molarRatios = c(N=1)),
    list(compoundName = "NO2", molarRatios = c(N=1)),
    list(compoundName = "N2O", molarRatios = c(N=1)),
    list(compoundName = "O2", molarRatios = c(O=1)),
    list(compoundName = "SO4", molarRatios = c(S=1)),
    list(compoundName = "CO2", molarRatios = c(C=1), sourceSink = T),
    list(compoundName = "N2", molarRatios = c(N=1), sourceSink = T),
    list(compoundName = "X", molarRatios = c(O=1), sourceSink = T),
    list(compoundName = "H2S", molarRatios = c(S=1), sourceSink = T)
  )

  gangstas = lapply(compoundParams, do.call, what = compoundFactory)

  processParams = list(
    list(
      processName = "AssimDOM",
      energyTerm = -4.32E-04,
      fromCompoundNames = list(C = "DOM", N = "DOM"),
      toCompoundNames = list(C = "Het", N = "Het"),
      massTerm = list(C = 1, N = 16/106),
      organismNames = c("Het")
    ),
    list(
      processName = "AssimCO2",
      energyTerm = -3.5E-02,
      fromCompoundNames = list(C = "CO2"),
      toCompoundNames = list(C = "Aut"),
      massTerm = list(C = 1),
      organismNames = c("Aut")
    ),
    list(
      processName = "AssimCH4",
      energyTerm = -1.09E-03,
      fromCompoundNames = list(C = "CH4"),
      toCompoundNames = list(C = "Met"),
      massTerm = list(C = 1),
      organismNames = c("Met")
    ),
    list(
      processName = "AssimNO3",
      energyTerm = -1.55E-04,
      fromCompoundNames = list(N = "NO3"),
      toCompoundNames = list(N = c("Het", "Aut", "Met")),
      massTerm = list(C = 1),
      organismNames = c("Het", "Aut", "Met")
    )
  )

##### PROBLEM: AssimNO3 will create 9 transformations: NO3 to Het, Aut, and Met for EACH ORGANISM!
  ## solution?  Rather than repeat for each organism, maybe match/recycle organism with from, to, mass.
  ## This would mean entering organismNames = c("Het", "Het") for catabolic processes like Denitrification
  ## where carbon source can be DOM or Biomass.



  return(unlist(gangstas, recursive = F))
}

#processFactory = function(gangstaObjects, processName, energyTerm, fromCompoundNames, toCompoundNames, massTerms, organismNames = "") {


processParams = list(
  list(
    processName = "AssimNO3",
    attributeList =
      list(
        energy = -1.55E-04,
        transformationIDs = c("toHet", "toAut", "toMet"),
        transformations = list(
          N = list(from = "NO3", to = c("Hetrph", "Autrph", "Metrph"), umols = 1.0)
        )
      )
  ),
  list(
    processName = "AssimNO2",
    attributeList =
      list(
        energy = -1.25E-04,
        transformationIDs = c("toHet", "toAut", "toMet"),
        transformations = list(
          N = list(from = "NO2", to = c("Hetrph", "Autrph", "Metrph"), umols = 1.0)
        )
      )
  ),
  list(
    processName = "AssimNH4",
    attributeList =
      list(
        energy = -3.18E-05,
        transformationIDs = c("toHet", "toAut", "toMet"),
        transformations = list(
          N = list(from = "NH4", to = c("Hetrph", "Autrph", "Metrph"), umols = 1.0)
        )
      )
  ),
  list(
    processName = "FeReduction",
    attributeList =
      list(
        organism = "Hetrph",
        energy = 3.0E-05,
        transformationIDs = c("ofDOM", "ofHetC"),
        transformations = list(
          C = list(from = c("DOM", "Hetrph"), to = NA, umols = 1.0),
          Fe = list(from = c("Fe3", "Fe3"), to = "Fe2", umols = 2.0),
          N = list(from = c("DOM", "Hetrph"), to = "NH4", umols = NA),
          Fe = list(from = c("DOM", "Hetrph"), to = "Fe2", umols = NA)
        )
      )
  ),
  list(
    processName = "Aerobic",
    attributeList =
      list(
        organism = "Hetrph",
        energy = 4.37E-04,
        transformationIDs = c("ofDOM", "ofHetC"),
        transformations = list(
          C = list(from = c("DOM", "Hetrph"),  to = NA,    umols = 1.0),
          O = list(from = c("O2",  "O2"    ),  to = NA,    umols = 2.0),
          N = list(from = c("DOM", "Hetrph"),  to = "NH4", umols = NA),
          Fe = list(from = c("DOM", "Hetrph"), to = "Fe2", umols = NA)
        )
      )
  ),
  list(
    processName = "DenitNO3toNO2",
    attributeList = list(
      organism = "Hetrph",
      energy = 2.88E-04,
      transformationIDs = c("ofDOM", "ofHetC"),
      transformations = list(
        C = list(from = c("DOM", "Hetrph"),  to = NA,   umols = 1.0),
        N = list(from = c("NO3", "NO3"   ),  to = "NO2",  umols = 2.0),
        N = list(from = c("DOM", "Hetrph"),  to = "NH4", umols = NA),
        Fe = list(from = c("DOM", "Hetrph"), to = "Fe2", umols = NA)
      )
    )
  ),
  list(
    processName = "DenitNO2toN2O",
    attributeList = list(
      organism = "Hetrph",
      energy = 4.15E-04,
      transformationIDs = c("ofDOM", "ofHetC"),
      transformations = list(
        C = list(from = c("DOM", "Hetrph"),  to = NA,   umols = 1.0),
        N = list(from = c("NO2", "NO2"   ),  to = "N2O",  umols = 2.0),
        N = list(from = c("DOM", "Hetrph"),  to = "NH4", umols = NA),
        Fe = list(from = c("DOM", "Hetrph"), to = "Fe2", umols = NA)
      )
    )
  ),
  list(
    processName = "DenitN2OtoN2",
    attributeList = list(
      organism = "Hetrph",
      energy = 6.45E-04,
      transformationIDs = c("ofDOM", "ofHetC"),
      transformations = list(
        C = list(from = c("DOM", "Hetrph"),  to = NA,   umols = 1.0),
        N = list(from = c("N2O", "N2O"   ),  to = NA,   umols = 4.0),
        N = list(from = c("DOM", "Hetrph"),  to = "NH4", umols = NA),
        Fe = list(from = c("DOM", "Hetrph"), to = "Fe2", umols = NA)
      )
    )
  ),
  list(
    processName = "SulfateRed",
    attributeList = list(
      organism = "Hetrph",
      energy = 3.8E-05,
      transformationIDs = c("ofDOM", "ofHetC"),
      transformations = list(
        C = list(from = c("DOM", "Hetrph"),  to = NA,    umols = 1.0),
        S = list(from = c("SO4", "SO4"   ),  to = NA,    umols = 0.5),
        N = list(from = c("DOM", "Hetrph"),  to = "NH4", umols = NA),
        Fe = list(from = c("DOM", "Hetrph"), to = "Fe2", umols = NA)
      )
    )
  ),
  list(
    processName = "Methanogen",
    attributeList = list(
      organism = "Hetrph",
      energy = 2.8E-05,
      transformationIDs = c("ofDOM", "ofHetC"),
      transformations = list(
        C = list(from = c("DOM", "Hetrph"),  to = NA,    umols = 0.5),
        C = list(from = c("DOM", "Hetrph"),  to = "CH4",   umols = 0.5),
        N = list(from = c("DOM", "Hetrph"),  to = "NH4", umols = NA),
        Fe = list(from = c("DOM", "Hetrph"), to = "Fe2", umols = NA)
      )
    )
  ),
  list(
    processName = "NitrifNH4toNO2",
    attributeList = list(
      organism = "Autrph",
      energy = 1.83E-04,
      transformations = list(
        N = list(from = "NH4",  to = "NO2",    umols = 2/3),
        O = list(from = "O2",   to = NA,     umols = 2.0)
      )
    )
  ),
  list(
    processName = "NitrifNO2toNO3",
    attributeList = list(
      organism = "Autrph",
      energy = 1.48E-04,
      transformations = list(
        N = list(from = "NO2",  to = "NO3",    umols = 2.0),
        O = list(from = "O2",   to = NA,     umols = 2.0)
      )
    )
  ),
  list(
    processName = "MethaneOxid",
    attributeList = list(
      organism = "Metrph",
      energy = 4.09E-04,
      transformations = list(
        C = list(from = "CH4",  to = NA,       umols = 0.5),
        O = list(from = "O2" ,  to = NA,       umols = 2.0)
      )
    )
  )
)

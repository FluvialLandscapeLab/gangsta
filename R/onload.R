.onLoad <- function(libname, pkgname) {

  op <- options()
  op.gangsta <- list(
    gangsta.path = "~/R-dev",
    gangsta.install.args = "",
    gangsta.name = "Geoffrey Poole",
    gangsta.desc.author = '"Geoffrey Poole <gpoole@montana.edu> [aut, cre]"',
    gangsta.desc.license = "What license is it under?",
    gangsta.desc.suggests = NULL,
    gangsta.desc = list(),

    gangsta.classes =
      c(
        base = "gangsta",
        comp = "compound",
        org = "organism",
        pool = "pool",
        proc = "process",
        metab = "metabolic",
        trans = "transformation"
      ),

    gangsta.vars =
      c(
        respEnergy = "respirationEnergyJoules",
        respRate = "RespRate",

        startSuffix = "initialAmountMols",
        endSuffix = "finalAmountMols",
        energySuffix = "energyJoules",
        transSuffix = "amountMolsTransfer"
      ),

    gangsta.attributes = c(
      name = "name",
      respRate = "respirationRate",
      orgName = "organismName",
      procName = "processName",
      compName = "compoundName",
      joulesToMols = "joulesToMolsRatio",
      fromPool = "from",
      toPool = "to",
      limitToStartMols = "limitToInitMols",
      molRatio = "molarRatio",
      initMols = "initialMols",
      finalMols = "finalMols",
      sourceSink = "sourceSink",
      energy = "energyTerm",
      element = "elementName",
      transOptions = "transferOptions"
    )


  )
#  toset <- !(names(op.gangsta) %in% names(op))
#  if(any(toset)) options(op.gangsta[toset])
   options(op.gangsta)
  invisible()
}



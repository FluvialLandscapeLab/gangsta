.onLoad <- function(libname, pkgname) {

  op <- options()
  op.gangsta <- list(
    gangsta.path = "~/R-dev",
    gangsta.install.args = "",
    gangsta.name = "Geoffrey Poole", ## What about Ann Marie?  She did some work too ;-)
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
        respEnergy = "respirationEnergy",
        respRate = "respirationRate",

        startSuffixPool = "initialAtoms",
        endSuffixPool = "finalAtoms",

        startSuffixCompound = "initialMolecules",
        endSuffixCompound = "finalMolecules",

        energySuffixProcess = "netEnergy",
        energySuffixOrganism = "totalProcessEnergy",

        transSuffix = "atoms"
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



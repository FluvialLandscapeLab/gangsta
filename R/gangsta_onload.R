.onLoad <- function(libname, pkgname) {

  op <- options()
  op.gangsta <- list(
    gangsta.path = "~/R-dev",
    gangsta.install.args = "",
    gangsta.name = c("Geoffrey Poole", "Ann Marie Reinhold"),
    gangsta.desc.author =
      c(
        '"Geoffrey Poole <gpoole@montana.edu> [aut, cre]"',
        '"Ann Marie Reinhold <annmarie.reinhold@montana.edu> [aut, cre]"'
      ),
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
        trans = "transfer"
      ),

    gangsta.vars =
      c(
        respEnergy = "respirationEnergy",
        respRate = "respirationRate",

        startSuffixPool = "initialAtoms",
        startSuffixCompound = "initialMolecules",
        startSuffixIsotopicRatio = "initialIsotopicRatio",

        endSuffixPool = "finalAtoms",
        endSuffixCompound = "finalMolecules",
        endSuffixIsotopicRatio = "finalIsotopicMols",

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

      molarAffinity = "molarAffinity",
      fromPool = "from",
      toPool = "to",
      limitToStartMols = "limitToInitMolecules",
      molRatio = "molarRatio",
      initialMolecules = "initialMolecules",
      finalMolecules = "finalMolecules",
      isotopicRatios = "isotopicRatios",
      infiniteCompound = "infiniteCompound",
      energy = "energyTerm",
      element = "elementName",
      transOptions = "transferOptions"
    )
  )
  options(op.gangsta)
  invisible()
}



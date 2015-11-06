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
        bnd = "bound",
        proc = "process",
        metab = "metabolic",
        trans = "transformation"
      ),

    gangsta.vars =
      c(
        respEnergy = "RespEnergy",
        respRate = "RespRate",

        startSuffix = "InitialMass",
        endSuffix = "FinalMass",
        energySuffix = "Energy",
        massSuffix = "Mass",
        transSuffix = "MassTrans"
      ),

    gangsta.attributes = c(
      name = "name",
      respRate = "respirationRate",
      refPool = "referencePoolName",
      orgName = "organismName",
      procName = "processName",
      compName = "compoundName",
      joulesToMols = "joulesToMolsRatio",
      fromPool = "from",
      toPool = "to",
      limitToStartMass = "limitToInitMols",
      molRatio = "molarRatio",
      initMols = "initialMols",
      sourceSink = "sourceSink",
      energy = "energyTerm",
      element = "elementName"
    )


  )
#  toset <- !(names(op.gangsta) %in% names(op))
#  if(any(toset)) options(op.gangsta[toset])
   options(op.gangsta)
  invisible()
}



.onLoad <- function(libname, pkgname) {

  op <- options()
  op.gangsta <- list(
    gangsta.path = "~/R-dev",
    gangsta.install.args = "",
    gangsta.name = "Geoffrey Poole",  ## what about AM Reinhold?
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
        metab = "metabolic", ## Metabolic is not a class in the current version of the model.
        trans = "transformation"
      ),

    gangsta.vars =  ## according to the new UML, all of these are attributes
      c(
        respEnergy = "respirationEnergy",
        respRate = "respirationRate",

        startSuffixPool = "initialAtoms",
        startSuffixCompound = "initialMolecules",

        endSuffixPool = "finalAtoms",
        endSuffixCompound = "finalMolecules",

        energySuffixProcess = "netEnergy",
        energySuffixOrganism = "totalProcessEnergy",

        transSuffix = "atoms"
      ),

    gangsta.attributes = c(
      name = "name",
      respRate = "respirationRate",

      orgName = "organismName",  ## should these next 3 just be "name" or "Organism.name", "Process.name", and "Compound.name"
      procName = "processName",
      compName = "compoundName",

      molarAffinity = "molarAffinity",
      fromPool = "from",
      toPool = "to",
      limitToStartMols = "limitToInitMolecules",  ## I think we need to add this to our UML diagram
      molRatio = "molarRatio",
      initialMolecules = "initialMolecules",
      finalMolecules = "finalMolecules",
      InfiniteCompound = "InfiniteCompound",  ## This is actually a class
      energy = "energyTerm",
      element = "elementName", ## I think we need to add this to our UML diagram
      transOptions = "transferOptions"  ## I think that we should change this to "multiToPools"
    )


  )
#  toset <- !(names(op.gangsta) %in% names(op))
#  if(any(toset)) options(op.gangsta[toset])
   options(op.gangsta)
  invisible()
}



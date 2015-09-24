#' Create GANGSTA Objects
#'
#' \code{compoundFactory} and \code{processFactory} are the primary functions
#' used to create GANGSTA objects.  \code{compoundFactory} creates both compound
#' and pool objects and \code{processFactory} creates processes and
#' transformation objects.  The constructors \code{compound}, \code{pool},
#' \code{process}, and \code{transformation} can be called individually, but
#' this is discouraged since the factory functions do substantial error checking
#' and assure that references between the GANGSTA objects are correct.
#'
#' \code{compoundFactory} is a function that creates a \code{compound} object
#' and all associated \code{pool} objects for use by the GANGSTA system.  It
#' calls the constructor methods for \code{compound} and \code{pool} objects.
#'
#' \code{processFactory} is similarly the prefered way to create \code{process}
#' objects and all assocaited \code{transformation} objects by calling the
#' constructor methods for \code{process} and \code{transformation} objects
#'
#' Direct use of constructors for \code{pool}, \code{compound}, \code{process},
#' and \code{transformation} should not be necessary and is discouraged.  Access
#' to these constructors is provided only for extensibility.
#'
#' GANGSTA uses several classes of S3 objects to represent biogeochemical
#' systems.  All S3 objects in GANGSTA built atop named lists with the class
#' arribute set.  The list's names are used as attribute names and the values in
#' the lists are the attribute values.  Thus, attribute values of the GANGSTA S3
#' objects are accessible with the notation x$name.
#'
#' \code{Compound} objects represent chemical species (e.g., SO4 and H2S for
#' sulfur) that contain one or more chemical elements; the mols (or umols, etc.)
#' of each element of interest within a \code{compound} are tracked using a
#' \code{pool}.  Each \code{pool} records the mass of an element in the
#' \code{compound} at the beginning and end a simulation timestep.  One
#' \code{pool} in the \code{compound} is designated as the "reference
#' \code{pool}". The molarRatios of all other \code{pools} in the compound are
#' relative to the reference \code{pool}.  The element of the reference
#' \code{pool} must appear first in the \code{molarRatios} vector, and the first
#' value in said vector must be 1.0 (i.e., the molar ratio of the reference
#' element to itself).
#'
#' Not all elements in a \code{compound} must be tracked.  The model developer
#' only creates \code{pool}s for the elements of interest.  For instance, when
#' creating a \code{compound} for SO4, passing \code{c(S=1)} to
#' \code{compoundFactory} as the \code{molarRatios} vector would create only a
#' \code{pool} for sulfur, not for oxygen.  Two \code{pool}s, one for sulfur and
#' one for oxygen, would be created by passing \code{c(S=1, O=4)}.  In both
#' cases, sulfur would be the reference \code{pool}.  To make oxygen the
#' reference \code{pool}, \code{molarRatios} would be \code{c(O = 1.0, S =
#' 0.25)}.
#'
#' \code{Compound}s have attributes named \code{name} and
#' \code{referencePoolName}. \code{Compound} objects are also used to represent
#' organisms (which assimilate elements as they grow) in GANGSTA models.
#' \code{Organism} objects inherit from \code{compound} and contain an extra
#' attribute called \code{respirationRate} For \code{organisms}, the
#' \code{respirationRate} is per mol (or umol, etc.) of element in the reference
#' \code{pool} per unit time in the simulation.
#'
#' \code{Pool}s have attributes called \code{name}, \code{elementName}, and
#' \code{compoundName}.  \code{Pool}s that are not designated as "reference
#' \code{pool}s" (i.e., those \code{pool}s that track elements appearing second
#' or later in the \code{molarRatios} vector) are consider to be \code{bound
#' pool}s.  \code{Bound pool} objects inhert from \code{pool} objects, and have
#' an additional attribute called \code{molarRatio}.
#'
#' GANGSTA models can operate using any unit of atomic count unit (mols, umols,
#' etc.), unit of energy (Joules, KJ, etc.) over any time unit defined by the
#' user.  However, it is critical that all units for values passed to the model
#' be consistent.  Here, the units of \code{respirationRate} and the units of
#' atomic count used by \code{pool} objects must be consistent with the units of
#' all other values passed to functions in the GANGSTA package.
#'
#' @param compoundName Name of the \code{compound} to be created (or for
#'   \code{pool}, the name of the \code{compound} to which the \code{pool}
#'   belongs).
#' @param molarRatios A named vector.  Names are the names of the chemical
#'   elements (think 'periodic table in chemistry') that are in the
#'   \code{compound} and that are to be tracked in the GANGSTA model.  Values in
#'   the vector are the ratios of each element in the \code{compound} relative
#'   to the first value in the vector (the "reference element"); thus the first
#'   value in the vector must always be 1.0 (the molar ratio of the first
#'   element to itself).
#' @param respirationRate The respiration rate (in units of energy per atomic
#'   count per time-1).  Atomic count refers to the mols (or umols, etc.) of the
#'   reference element. When \code{respirationRate} is numeric, an
#'   \code{organism} object is returnd.  When NA, a \code{compound} object is
#'   returned.
#' @param sourceSink Boolean when set to TRUE tags a compound as being unlimited
#'   in supply for the purposes of the model.
#' @param gangstaObjects A list of compounds and pools, typically created by
#'   calling \code{compoundFactory}.  Error checking in \code{processFactory}
#'   checks to be sure that referenced compounds and pools exist in the
#'   gangstaObjects list.
#' @param processName. The name of the process to be created.
#' @param energyTerm. The net energy associated with the processs.  A positive
#'   number represents a process that yeilds energy, a negative number
#'   represents a process that consumes energy.  Units are kJ (or J, etc.) of
#'   energy per mol (or umol, etc.) listed in \code{molarTerms} parameter.
#' @param fromCompoundNames Named \code{list} of compound names where the name
#'   of each \code{list} member is the a chemical element derived from the
#'   compound.  For instance, to track carbon flux from the oxidation of
#'   glucose, the \code{fromCompoundNames} list might be list(C = "C6H12O6", O =
#'   "O2")
#' @param toCompoundNames.  See \code{toCompoundNames}.  To track carbon flux
#'   from the oxidation of glucose, the \code{toCompoundNames} list might be
#'   list(C = "CO2", O = "Ox") (where Ox is a undifferentiated sink for oxygen
#'   comprised of H2O and CO2).  The names of \code{toCompoundNames} must be the
#'   same and in the same order as those of \code{fromCompoundNames}.
#' @param molarTerms Named list containing the mols (or uMols, etc.) of each
#'   element that are transformed by the process.  The names of \code{molarTerms}
#'   must be the same and in the same order as those of \code{fromCompoundNames}
#'   and \code{toCompoundNames}.
#' @param organismNames A vector of organisms that utilize the process.
#' @param elementName The name of the element contained by the created
#'   \code{pool}.
#' @param molarRatio The ratio of the elemental mass in a \code{bound pool} to
#'   the mass in its reference \code{pool}.  When molarRatio is NA, a
#'   \code{pool} object is return.  When molarRatio is numeric, a \code{bound
#'   pool} object of returned.
#' @param referencePoolName Name of the reference \code{pool} for the
#'   \code{compound} object to be created.
#' @return \code{compoundFactory} returns a list of \code{compound} and {pool}
#'   objects. \code{processFactory} return a list of \code{process} and
#'   \code{transformation} objects.  The remaining constructor methods return an
#'   individual GANGSTA object of the class corresponding to the function name.

compoundFactory = function(compoundName, molarRatios, initialMols, respirationRate = NA, sourceSink = F) {
  checkNames = unique(names(molarRatios))==""
  if(any(checkNames) || (length(checkNames) != length(molarRatios))) {
    stop("Each member of the molarRatios vector must be named using an element name.  Element names must be unique.")
  }
  elementNames = names(molarRatios)
  if(molarRatios[1]!=1.0) {
    stop("The molarRatio of the Reference Element (the first element in the molarRatios vector) must be 1.0")
  }
  molarRatios[1] = NA
  newPools = mapply(pool, compoundName, elementNames, molarRatios, USE.NAMES = F, SIMPLIFY = F)
  names(newPools) = sapply(newPools, function(x) x$name)
  newCompound = list(compound(compoundName, newPools[[1]]$name, initialMols, respirationRate, sourceSink))
  names(newCompound) = compoundName
  return(c(newCompound, newPools))
}

#' @rdname compoundFactory
processFactory = function(gangstaObjects, name, energyTerm, fromCompoundNames, toCompoundNames, molarTerms, organismName = "", limitToInitMols = T) {
  if(!identical(organismName, "")) {
    gangstasExist(gangstaObjects, organismName, "organism")
  }

  inputList = list(name, energyTerm, organismName, limitToInitMols)
  if(any(plyr::laply(inputList, length) != 1)) {
    stop(paste0("Process: ", name, "\n The arguments name, energyTerm, orgaismName, and limitToInitMols must be vectors of length() = 1."))
  }

  inputList = list(fromCompoundNames, toCompoundNames, molarTerms)
  if(length(unique(plyr::laply(inputList, length)))!=1) {
    stop(paste0("Process: ", name, "\n The length of fromCompoundNames, toCompoundNames, and molarTerms vectors must be equal."))
  }

  nullNames = plyr::laply(inputList, function(x) is.null(names(x)), .drop = F)
  if(any(nullNames)) {
    stop(stop(paste0("Process: ", name, "\n The members of lists fromCompoundNames, toCompoundNames, and molarTerms must be named.")))
  }

  elementMatrix = plyr::laply(inputList, names, .drop = F)

  differentNamesAcrossLists = apply(elementMatrix, 2, function(x) length(unique(x)) != 1)
  if(any(differentNamesAcrossLists)) {
    stop(paste0("Process: ", name,": \n The members of 'fromCompoundNames,' 'toCompoundNames,' and 'molarTerms' lists must have the same names in the same order across lists."))
  }

  newProcess = list(process(name, energyTerm, organismName))
  names(newProcess) = name

  fromCompoundNames = replaceDotWithOrganism(fromCompoundNames, organismName)
  toCompoundNames = replaceDotWithOrganism(toCompoundNames, organismName)
  fromPoolNames = makePoolNames(fromCompoundNames, gangstaObjects = gangstaObjects)
  toPoolNames = makePoolNames(toCompoundNames, gangstaObjects = gangstaObjects)
  molarTerms = replaceNAWithMolarRatio(molarTerms, fromPoolNames, fromCompoundNames, gangstaObjects)

  newTransformations = mapply(transformation, fromPoolNames, toPoolNames, molarTerms,
                              MoreArgs = list(gangstaObjects = c(gangstaObjects, newProcess), processName = name, limitToInitMols = limitToInitMols),
                              SIMPLIFY = F)
  names(newTransformations) = sapply(newTransformations, function(x) x$name)

  duplicateNames = unique(names(newTransformations)[duplicated(names(newTransformations))])
  if(length(duplicateNames)>0) {
    stop("The following transformation(s) were specified more than once: ", paste0(duplicateNames, collapse = "; "))
  }

  return(c(newProcess, newTransformations))
}

#' @rdname compoundFactory
compound = function(compoundName, referencePoolName, initialMols, respirationRate = NA, sourceSink) {
  newCompound = list(name = compoundName, referencePoolName = referencePoolName, initialMols = initialMols, sourceSink = sourceSink)
  class(newCompound) = c("compound", "gangsta")
  if(!is.na(respirationRate)) {
    newCompound = structure(c(newCompound, list(respirationRate = respirationRate)), class = c("organism", class(newCompound)))
  }
  return(newCompound)
}

#' @rdname compoundFactory
pool = function(compoundName, elementName, molarRatio = NA) {
  poolName = makePoolNames(compoundName, elementName)
  newPool = list(name = poolName, elementName = elementName, compoundName = compoundName)
  class(newPool) = c("pool", "gangsta")
  if(!is.na(molarRatio)) {
    newPool = structure(c(newPool, list(molarRatio = molarRatio)), class = c("bound", class(newPool)))
  }
  return(newPool)
}

#' @rdname compoundFactory
process = function(name, energyTerm, organismName = "") {
  processClassNames = c(gangstaClassName("proc"), gangstaClassName("base"))
  newProcess = list(name = name, energyTerm = energyTerm)
  class(newProcess) = processClassNames
  if(energyTerm != 0) {
    newProcess = structure(c(newProcess, list(organismName = organismName)), class = c(gangstaClassName("metab"), class(newProcess)))
  }
  return(newProcess)
}

#' @rdname compoundFactory
transformation = function(gangstaObjects, processName, fromPoolName, toPoolName, molarTerm, limitToInitMols = T){
  # Calling fromToPair does some key error checking.
  pools = fromToPair(gangstaObjects, fromPoolName, toPoolName)
  transformationName = paste(processName, fromPoolName, toPoolName, sep="_")
  process = getGangstas(gangstaObjects, processName)
  energyToMolsRatio = process[[1]]$energyTerm / molarTerm
  newTransformation =
    list(
      name = transformationName,
      from = fromPoolName,
      to = toPoolName,
      molarTerm = molarTerm,
      energyToMolsRatio = energyToMolsRatio,
      processName = processName,
      limitToInitMols = limitToInitMols
    )
  class(newTransformation) = c("transformation", "gangsta")
  return(newTransformation)
}

processSpec = function(listOfProcessFactoryArgs) {
  checkProcessSpecNames(names(listOfProcessFactoryArgs))
  class(listOfProcessFactoryArgs) = c("processSpec", "gangsta")
  return(listOfProcessFactoryArgs)
}


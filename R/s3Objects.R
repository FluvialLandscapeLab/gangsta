#' Create GANGSTA compounds and pools
#'
#' \code{compoundFactory} is a function that creates a \code{compound} object
#' and all associated \code{pool} objects for use by the GANGSTA system.  It
#' calls the constructor methods for \code{compound} and \code{pool} objects.
#' Direct use of constructors for \code{pool} and \code{compound} should not be
#' necessary.  Access to these constructors is provided only for extensibility.
#' Use \code{compoundFactory} to create \code{pool} and \code{compound} objects.
#'
#' GANGSTA uses several classes of S3 objects to represent biogeochemical
#' systems.  All S3 objects in GANGSTA are named lists with the class arribute
#' set.  The list's names are used as attribute names and the values in the
#' lists are the attribute values.  Thus, attribute values of the GANGSTA S3
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
#' @param compoundName Name of the \code{compound} to be created (or for \code{pool},
#'   the name of the \code{compound} to which the \code{pool} belongs).
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
#' @param elementName The name of the element contained by the created
#'   \code{pool}.
#' @param molarRatio The ratio of the elemental mass in a \code{bound pool} to
#'   the mass in its reference \code{pool}.  When molarRatio is NA, a
#'   \code{pool} object is return.  When molarRatio is numeric, a \code{bound
#'   pool} object of returned.
#' @param referencePoolName Name of the reference \code{pool} for the
#'   \code{compound} object to be created.
#' @return \code{compoundFactory} returns a list of length 2.  The first item in
#'   the list is the new \code{compound} (or \code{organism}) object.  The
#'   second is a list of the associated \code{pool} (and \code{bound pool})
#'   objects.

compoundFactory = function(compoundName, molarRatios, respirationRate = NA, sourceSink = F) {
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
  newCompound = list(compound(compoundName, newPools[[1]]$name, respirationRate, sourceSink))
  names(newCompound) = compoundName
  return(c(newCompound, newPools))
}

#' @rdname compoundFactory
pool = function(compoundName, elementName, molarRatio = NA) {
  poolName = makePoolName(compoundName, elementName)
  newPool = list(name = poolName, elementName = elementName, compoundName = compoundName)
  class(newPool) = c("pool", "gangsta")
  if(!is.na(molarRatio)) {
    newPool = structure(c(newPool, list(molarRatio = molarRatio)), class = c("bound", class(newPool)))
  }
  return(newPool)
}

#' @rdname compoundFactory
compound = function(compoundName, referencePoolName, respirationRate = NA, sourceSink) {
  newCompound = list(name = compoundName, referencePoolName = referencePoolName, sourceSink = sourceSink)
  class(newCompound) = c("compound", "gangsta")
  if(!is.na(respirationRate)) {
    newCompound = structure(c(newCompound, list(respirationRate = respirationRate)), class = c("organism", class(newCompound)))
  }
  return(newCompound)
}



# HetAerobicResp =
#   list(
#     fromCompounds = list(C = c("DOM", "Het"), O = c("O2")),
#     toCompounds = list(C = "CO2", O = "X"),
#     massTerms = list(C = 1, O = 2),
#     organisms = list("Het")
#   )
#

#      fromCompounds = list(N = c("DOM", "Het"))
#      toCompounds = list(N = "NH4")
#      massTerms = list(N = 1)
#      organisms = list("Het")

processFactory = function(gangstaObjects, processName, energyTerm, fromCompoundNames, toCompoundNames, massTerms, organismNames = "") {

  if(!identical(organismNames, "")) {
    gangstasExist(gangstaObjects, organismNames, "organism")
  }

  badNames = F
  elementMatrix = sapply(list(fromCompoundNames, toCompoundNames, massTerms), function(inputList) suppressWarnings(unique(names(inputList))))
  ## sapply above will return a matrix if all the lists are the same length.
  if(!is.matrix(elementMatrix)){
    badNames = T
  } else {
    ## if names are the same,
    nonUniqueElements = apply(elementMatrix, 1, function(x) length(unique(x)) != 1)
    if(any(nonUniqueElements)) {
      badNames = T
    }
  }

  if(badNames) {
    stop("Parameters 'fromCompoundNames,' 'toCompoundNames,' and 'massTerms' must be named vectors. \n  - Names must be chemical elements, which must be the same and in the same order across lists, but unique within each list.\n  - All list  the same for each list.\n  - No chemical ename can be duplicated within a vector.")
  }

  processNames = paste0(organismNames, processName)
  newProcesses = mapply(process, processName = processNames, energyTerm = energyTerm, organismName = organismNames, SIMPLIFY = F)

  fromPoolNames = makeMultiplePoolNames(fromCompoundNames)
  toPoolNames = makeMultiplePoolNames(toCompoundNames)

  gangstasExist(gangstaObjects, unlist(fromPoolNames), "pool")
  gangstasExist(gangstaObjects, unlist(toPoolNames), "pool")


  ## To understand this triple nested mapply, think of each mapply as a
  ## nested "for loop."
  ##
  ## The outermost mapply repeats for each process in "processNames" (which is
  ## always 1:1 with length(organismNames)).
  ##
  ## For each process, then, the next mapply repeats for each chemical element
  ## in the fromPools, toPools, and massTerms lists.
  ##
  ## For each process and element, then, the innermost loop matches (with
  ## recycling) the toPools, fromPools, and massTerms and calls the
  ## transformation() function to make a transformation for each fromPool,
  ## toPool, and massTerm triplet.
  ##
  ## The strategically placed "unlist" functions yield a flat list of the
  ## resulting transformations, rather than transformations clustered into
  ## sublists by process and chemical element.
  newTransformations =
##    unlist(
##      mapply(
##        function(process, org)
          unlist(
            mapply(
              function(froms, tos, mTerms)
                mapply(
                  transformation,
                  froms,
                  tos,
                  mTerms,
                  MoreArgs = list(
                    processName = processNames,
                    gangstaObjects = c(gangstaObjects, newProcesses)
                  ),
                  SIMPLIFY = F
                ),
              fromPoolNames,
              toPoolNames,
              massTerms
            ),
            recursive = F
##          ),
##        processNames,
##        organismNames,
##        SIMPLIFY = F
##      ),
##      recursive = F
    )

  return(c(newProcesses, newTransformations))
}

process = function(processName, energyTerm, organismName = NA) {
  if(is.na(organismName)) {
    organismName = ""
  }
  newProcess = list(name = processName, energyTerm = energyTerm)
  class(newProcess) = c("process", "gangsta")
  if(organismName != "") {
    newProcess = structure(c(newProcess, list(organismName = organismName)), class = c("metabolic", class(newProcess)))
  }
  return(newProcess)
}

transformation = function(gangstaObjects, processName, fromPoolName, toPoolName, massTerm){
  # Calling fromToPair does some key error checking.
  pools = fromToPair(gangstaObjects, fromPoolName, toPoolName)
  transformationName = paste(processName, fromPoolName, toPoolName, sep="_")
  process = getGangstas(gangstaObjects, processName)
  energyToMassRatio = process[[1]]$energyTerm / massTerm
  newTransformation =
    list(
      name = transformationName,
      from = fromPoolName,
      to = toPoolName,
      massTerm = massTerm,
      energyToMassRatio = energyToMassRatio,
      processName = processName
    )
  class(newTransformation) = c("transformation", "gangsta")
  return(newTransformation)
}


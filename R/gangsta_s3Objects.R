# Attribute values of the \code{gangsta} S3 objects are accessible with
# the notation x$name. The constructors \code{compound}, \code{pool},
# \code{process}, and \code{transfer} can be called individually, but this is
# discouraged since the factory functions do substantial error checking and
# assure that references between the \code{gangsta} objects are correct.

#' Create \code{gangsta} objects
#'
#' \code{compoundFactory} and \code{processFactory} are the primary functions
#' used to create gangsta objects.  \code{compoundFactory} is a function that
#' creates a \code{compound} object and all associated \code{pool} objects for
#' use by the \code{gangsta} package. \code{processFactory} creates a
#' \code{process} object and all associated \code{transfer} objects.
#'
#' The \code{gangsta} uses several classes of S3 objects to represent
#' biogeochemical systems.  All S3 objects in \code{gangsta} are built atop
#' named lists with the class arribute set.  The list's names are used as
#' attribute names and the values in the lists are the attribute values.
#'
#' \code{compound} objects represent chemical species (e.g., SO4 and HS) that
#' generally contain one or more chemical elements; the mols (or umols, etc.) of
#' each element of interest (i.e., tracked element) within a \code{compound} are
#' tracked using a \code{pool}.  Each \code{pool} records the mols of a tracked
#' element in the \code{compound}.
#'
#' The \code{molarRatios} describe the chemical composition of the tracked
#' elements contained in compounds. For example, the \code{molarRatios} for
#' water could be represented as \code{c(H = 2, O = 1)} because a molecule of
#' water is composed of two hydrogen atoms and one oxygen atom.
#'
#' Not all elements in a \code{compound} must be tracked.  The model developer
#' only creates \code{pool}s for the elements of interest.  For instance, when
#' creating a \code{compound} for SO4, passing \code{c(S=1)} to
#' \code{compoundFactory} as the \code{molarRatios} vector would create only a
#' \code{pool} for sulfur, not for oxygen.  Two \code{pool}s, one for sulfur and
#' one for oxygen, would be created by passing \code{c(S=1, O=4)}.
#'
#' \code{compound} objects have attributes named \code{name},
#' \code{initialMolecules} which are the starting values for the model run, and
#' \code{infiniteCompound} which is a Boolean indicating whether the compound is
#' a source/sink.
#'
#' \code{compound} objects are also used to represent organisms (which
#' assimilate elements as they grow) in \code{gangsta}-derived models.
#' \code{organism} objects inherit from \code{compound} and contain an extra
#' attribute called \code{respirationRate} For \code{organism}s, the
#' \code{respirationRate} is energy (J, KJ, etc.) per mol (or umol, etc.) of
#' \code{compound} per timestep.
#'
#' \code{pool} objects have attributes called \code{name}, \code{elementName}
#' which is the chemical element stored in the \code{pool}, \code{compoundName}
#' which indicates the compound with which the pool is associated, and
#' \code{molarRatio} which describes the ratio of atoms in the \code{pool} per
#' unit molecule, i.e., \code{compound}.
#'
#' \code{process} objects have attributes called \code{name}, \code{energyTerm}
#' which is equal to the energy generated or consumed by the \code{process},
#' \code{organismName} which is the name of the \code{organism} carrying out the
#' \code{process}, and \code{transferOptions} which are indices for the
#' \code{transfer}s associated with each \code{process}.  Each \code{process}
#' with nonzero \code{energyTerm} is of class "metabolic."
#'
#' \code{transfer} objects have attributes called \code{name}, \code{from} which
#' is the \code{pool} from which atoms are being transfered, \code{to} which is
#' the \code{pool} to which atoms are being transfered, \code{molarTerm} which
#' enforces the stoichiometry of the \code{process}, \code{molarAffinity} which
#' is equal to the \code{energyTerm} of the \code{process} divided by the
#' \code{molarTerm} of the \code{transfer}, \code{processName} which is the name
#' of the \code{process} that the \code{transfer} is associated with, and
#' \code{limitToInitMolecules} which is a Boolean that inherits from
#' \code{process}.
#'
#' \code{gangsta}-derived models can operate using any unit of atomic count unit
#' (mols, umols, etc.), unit of energy (Joules, KJ, etc.) over any time unit
#' defined by the user.  However, it is critical that all units for values
#' passed to the model be consistent.  For example, the units of
#' \code{respirationRate} and the units of atomic count used by \code{pool}
#' objects must be consistent with the units of all other values passed to
#' functions in the \code{gangsta} package.
#'
#' @param compoundName A character vector of \code{length = 1} containing the
#'   name of the \code{compound} to be created.
#' @param molarRatios A named numeric vector.  Vector names are the names of the
#'   chemical elements (think 'periodic table in chemistry') that are in the
#'   \code{compound} and that are to be tracked by \code{gangsta}.  Values in
#'   the vector are the ratios for the number of atoms of each element in the
#'   \code{compound}.
#' @param initialMolecules The number of mols (or umols, etc.) of the compound
#'   available at the beginning of the simulation.  If the \code{compound} is of
#'   type \code{infiniteCompound}, then \code{initialMolecules} must be set to
#'   0.
#' @param respirationRate The respiration rate (in units of energy per unit of
#'   biomass per model timestep); applies to \code{organism} objects.  When
#'   \code{respirationRate} is numeric, an \code{organism} object is returned.
#'   When NA, a \code{compound} object is returned.  \code{respirationRate}
#'   values must be negative (i.e., energy cost to the organisms).
#' @param infiniteCompound Boolean when set to \code{TRUE} tags a
#'   \code{compound} as being unlimited in supply. \code{infiniteCompound}s
#'   represent sources/sinks.  When \code{infiniteCompound} is set to
#'   \code{TRUE}, \code{initialMolecules} must be set to 0.
#' @param gangstaObjects A list of objects of class \code{gangsta} representing
#'   the \code{organisms}, \code{compounds}, \code{pools}, \code{processes}, and
#'   \code{transfers} to be included in the model.  These objects are created
#'   using \code{compoundFactory} and \code{processFactory}.
#'   \code{compoundFactory} must be executed before \code{processFactory}
#'   because \code{gangstaObjects} for all \code{compounds} involved in a
#'   \code{process} must be created before \code{processFactory} can create the
#'   \code{process}.
#' @param processName A character vector of the name of the \code{process} to be
#'   created.
#' @param energyTerm The chemical affinity of the processs.  A positive number
#'   represents a process that yeilds energy, a negative number represents a
#'   process that consumes energy.  Units are kJ (or J, etc.) of energy per mol
#'   (or umol, etc.) of the reaction.
#' @param fromCompoundNames Named \code{list} where the name of each \code{list}
#'   member is the a chemical element derived from the compound.  For instance,
#'   to track carbon flux from the oxidation of glucose, the
#'   \code{fromCompoundNames} list might be \code{list(C = "C6H12O6", O = "O2")}
#' @param toCompoundNames  See \code{fromCompoundNames}.  To track carbon flux
#'   from the oxidation of glucose, the \code{toCompoundNames} list might be
#'   \code{list(C = "CO2", O = "Ox")} (where Ox is a undifferentiated sink for
#'   oxygen comprised of H2O and CO2).  The names of \code{toCompoundNames} must
#'   be the same and in the same order as those of \code{fromCompoundNames}.
#' @param molarTerms Named list containing the mols (or umols, etc.) of each
#'   element that are transfered by the \code{process}.  The names of
#'   \code{molarTerms} must be the same and in the same order as those of
#'   \code{fromCompoundNames} and \code{toCompoundNames}.
#' @param transferOptions When transferOptions is \code{NULL},
#'   \code{processFactory} will create these automatically.
#'   \code{transferOptions} consist of a list of integer or numeric vectors
#'   containing the indicies of transfers in a process.  This specification is
#'   appropriate when each \code{from pool} has exactly one \code{to pool}.
#'   However, these must be specified when more than one \code{to pool} exists
#'   for any transfer involved in a process, i.e., indicies are grouped when
#'   transfers represent optional pathways. For instance, if a process has four
#'   transfers (fromA -> toA, fromB -> toB1, fromB -> toB2, fromC -> toC), the
#'   second and third transfers can represent an option.  fromB can go to either
#'   toB1 or toB2, so long as the sum of the two options is in stoichiometric
#'   balance with the A and C tranfers.  To represent such an option, the
#'   transferOption list would be \code{list(1, 2:3, 4)}.
#' @param organismName Name of the organism carrying out the \code{process}.
#' @param elementList A named list of character vectors where each vector
#'   contains the names of the isotopes to be tracked for a single element.
#'   Typically, these are character representations of each isotope's molecular
#'   weight. The name of each element in this list corresponds to the element
#'   names used in the \code{molarRatios} list.
#'
#' @return \code{compoundFactory} returns a list of \code{compound} and
#'   \code{pool} objects. \code{processFactory} return a list of \code{process}
#'   and \code{transfer} objects.
#' @export
compoundFactory = function(compoundName, molarRatios, initialMolecules, respirationRate = NA, infiniteCompound = F) {
  checkNames = unique(names(molarRatios))==""
  if(any(checkNames) || (length(checkNames) != length(molarRatios))) {
    stop("Each member of the molarRatios vector must be named using an element name.  Element names must be unique.")
  }
  badInitialMolecules = infiniteCompound && (initialMolecules != 0)
  if(badInitialMolecules) stop("compoundFactory error for compound = ", compoundName, "; initialMolecules must be 0 for infiniteCompounds.")
  elementNames = names(molarRatios)
  newPools = mapply(pool, compoundName, elementNames, molarRatios, USE.NAMES = F, SIMPLIFY = F)
  names(newPools) = sapply(newPools, function(x) x$name)
  newCompound = list(compound(compoundName, initialMolecules, respirationRate, infiniteCompound))
  names(newCompound) = compoundName
  return(c(newCompound, newPools))
}

#' @rdname compoundFactory
#' @export
enableIsotopeTracking = function(gangstaObjects, elementList, initialIsotopicRatios){
  poolObjectIdx = subsetGangstas(gangstaObjects, "class", getOption("gangsta.classes")["pool"], asIndex = T)
  poolObjects = gangstaObjects[poolObjectIdx]
  elementNames = names(elementList)
  # error check for duplicate elements
  if(!(length(elementNames) == length(unique(elementNames)))){
    stop("One or more elementList names are duplicated")
  }
  # Sort pool objects to be replaced by element
  poolObjectsByElement = lapply(elementNames,
                                subsetGangstas,
                                gangstaObjects = poolObjects,
                                attributeName = "elementName")

  # Throw a warning for elements in element list that aren't in the model
  doesNotMatch = sapply(poolObjectsByElement, function(x) length(x)==0)
  if (any(doesNotMatch)) warning("The following elements are in elementList but not in the model:",
                              paste0(elementNames[doesNotMatch], collapse = ", " ))

  # Make new pools for each isotope
  isotopeSpecificPools =
    unlist(
      unlist(
        Map(function(poolObjects, isotopeNames){
          lapply(poolObjects,
                 function(poolObject){
                   isotopeSpecificPoolObjects = Map(pool,
                                                    compoundName = poolObject$compoundName,
                                                    elementName = paste0(poolObject$elementName,isotopeNames),
                                                    molarRatio = poolObject$molarRatio,
                                                    isotopicRatio = initialIsotopicRatios[[poolObject$name]])
                   return(isotopeSpecificPoolObjects)
                 }
          )
        },
        poolObjectsByElement,
        elementList
        ),
        recursive = F
      ),
      recursive = F
    )
  names(isotopeSpecificPools) = lapply(isotopeSpecificPools, function(poolObject) poolObject$name)

  # Throw an error if user did not include initialIsotopicRatios for a pool
  # associated with an element in ElementList

  # Throw an error if user supplied initialIsotopicRatios are not the
  # same length or don't have the same names as the isotopes for a particular
  # element specified in elementList

  # Make new list of gangsta objects by concatenating new isotope specific pools with
  # the old gangsta list subsetted to exclude pools that have been "split" by isotope

  newGangstaObjects = c(gangstaObjects[!names(gangstaObjects) %in% names(unlist(poolObjectsByElement, recursive = F))],
                        isotopeSpecificPools)

  return(newgangstaObjects)
}

#' @rdname compoundFactory
#' @export
processFactory = function(gangstaObjects, processName, energyTerm, fromCompoundNames, toCompoundNames, molarTerms, transferOptions = NULL, organismName = "", limitToInitMolecules = T) {

  # check to be sure the specifeid organism already exists in gangstaObjects
  if(!identical(organismName, "")) {
    gangstasExist(gangstaObjects, organismName, "organism")
  }

  # check to be sure some parameters are vectors of length 1
  inputList = list(processName, energyTerm, organismName, limitToInitMolecules)
  if(any(plyr::laply(inputList, length) != 1)) {
    stop(paste0("Process: ", processName, "\n The arguments processName, energyTerm, orgaismName, and limitToInitMolecules must be vectors of length() = 1."))
  }

  # check to be sure some vectors are equal in length
  inputList = list(fromCompoundNames, toCompoundNames, molarTerms)
  if(length(unique(plyr::laply(inputList, length)))!=1) {
    stop(paste0("Process: ", processName, "\n The length of fromCompoundNames, toCompoundNames, and molarTerms vectors must be equal."))
  }

  # check to be sure some required names are present
  nullNames = plyr::laply(inputList, function(x) is.null(names(x)), .drop = F)
  if(any(nullNames)) {
    stop(stop(paste0("Process: ", processName, "\n The members of lists fromCompoundNames, toCompoundNames, and molarTerms must be named.")))
  }

  # check to be sure from, to and molarTerms have same names
  elementMatrix = plyr::laply(inputList, names, .drop = F)
  differentNamesAcrossLists = apply(elementMatrix, 2, function(x) length(unique(x)) != 1)
  if(any(differentNamesAcrossLists)) {
    stop(paste0("Process: ", processName,": \n The members of 'fromCompoundNames,' 'toCompoundNames,' and 'molarTerms' lists must have the same names in the same order across lists."))
  }

  # create the process object
  if(is.null(transferOptions)) transferOptions = structure(as.list(1:length(fromCompoundNames)), names = names(fromCompoundNames))
  newProcess = list(process(processName, energyTerm, transferOptions, organismName))
  names(newProcess) = processName

  fromCompoundNames = replaceDotWithOrganism(fromCompoundNames, organismName)
  toCompoundNames = replaceDotWithOrganism(toCompoundNames, organismName)
  fromPoolNames = makePoolNames(fromCompoundNames, gangstaObjects = gangstaObjects)
  toPoolNames = makePoolNames(toCompoundNames, gangstaObjects = gangstaObjects)

  # expandedMolarTerms = rep.int(sapply(molarTerms, "[[", "value"), times = sapply(molarTerms, function(x) length(x$groupIdx)))
  # molarTerms = replaceNAWithMolarRatio(molarTerms, fromPoolNames, fromCompoundNames, gangstaObjects)

  newTransfers = mapply(transfer, fromPoolNames, toPoolNames, molarTerms,
                        MoreArgs = list(gangstaObjects = c(gangstaObjects, newProcess), processName = processName, limitToInitMolecules = limitToInitMolecules),
                        SIMPLIFY = F)
  names(newTransfers) = sapply(newTransfers, function(x) x$name)
  duplicateNames = unique(names(newTransfers)[duplicated(names(newTransfers))])
  if(length(duplicateNames)>0) {
    stop("The following transfer(s) were specified more than once: ", paste0(duplicateNames, collapse = "; "))
  }

  return(c(newProcess, newTransfers))
}

# @rdname compoundFactory
compound = function(compoundName, initialMolecules, respirationRate = NA, infiniteCompound) {
  newCompound = list(name = compoundName, initialMolecules = initialMolecules, infiniteCompound = infiniteCompound)
  class(newCompound) = c("compound", "gangsta")
  if(!is.na(respirationRate)) {
    if(respirationRate > 0) {
      stop("Respiration rate must be negative.")
    }
    newCompound = structure(c(newCompound, list(respirationRate = respirationRate)), class = c("organism", class(newCompound)))
  }
  return(newCompound)
}

# @rdname compoundFactory
# @param isotopeName The name of the isotope or element contained by the created
#   \code{pool}.
pool = function(compoundName, elementName, molarRatio, isotopicRatio = NA) {
  poolName = makePoolNames(compoundName, elementName)
  newPool = list(name = poolName, elementName = elementName, compoundName = compoundName, molarRatio = molarRatio)
  if(!is.na(isotopicRatio)){
    newPool$isotopicRatio = isotopicRatio
  }
  class(newPool) = c("pool", "gangsta")
  return(newPool)
}

# @rdname compoundFactory
element = function(elementName, isotopeNames){
  newElement = list(name = elementName,
                    isotopes = isotopeNames)
  class(newElement) = c("element", "gangsta")
  return(newElement)
}

# @rdname compoundFactory
process = function(processName, energyTerm, transferOptions, organismName = "") {
  processClassNames = c(gangstaClassName("proc"), gangstaClassName("base"))
  newProcess = list(name = processName, energyTerm = energyTerm, organismName = organismName, transferOptions = transferOptions)
  class(newProcess) = processClassNames
  if(energyTerm != 0) {
    class(newProcess) = c(gangstaClassName("metab"), class(newProcess))
  }
  return(newProcess)
}

# @rdname compoundFactory
transfer = function(gangstaObjects, processName, fromPoolName, toPoolName, molarTerm, limitToInitMolecules = T){
  # Calling fromToPair does some key error checking.
  pools = fromToPair(gangstaObjects, fromPoolName, toPoolName)
  transferName = paste(processName, fromPoolName, toPoolName, sep="_")
  process = getGangstas(gangstaObjects, processName)
  energyToMolsRatio = process[[1]]$energyTerm / molarTerm
  newTransfer =
    list(
      name = transferName,
      from = fromPoolName,
      to = toPoolName,
      molarTerm = molarTerm,
      molarAffinity = energyToMolsRatio,
      processName = processName,
      limitToInitMolecules = limitToInitMolecules
    )
  class(newTransfer) = c(gangstaClassName("trans"), "gangsta")
  return(newTransfer)
}




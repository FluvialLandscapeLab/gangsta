#' Query a list of \code{gangsta} objects
#'
#' Lists of \code{gangsta} objects created using \code{compoundFactory} and
#' \code{processFactory} can be queried conveniently using these methods.
#'
#' \code{subsetGangstas} Queries \code{gangsta} objects by attribute.  Returns a
#' subset of \code{gangstaObjects} where the attribute specified by
#' \code{attributeName} is equal to the value specified by
#' \code{attributeValue}.
#'
#' \code{gangstasExist} returns \code{TRUE} if all names in \code{gangstaNames}
#' have an associated \code{gangsta} object in \code{gangstaObjects}.  If not,
#' execution is halted and an error reported.
#'
#' \code{getGangstas} is a simplified wrapper for \code{subsetGangstas} and
#' queries \code{gangstaObjects} on the "name" attribute.  It also throws an
#' error if no object of a provided name exists in \code{gangstaObjects} or if
#' more than one object with a provided name exists in \code{gangstaObjects}.
#'
#' \code{getGangstaAttribute} returns a vector of attribute values from the
#' objects in \code{gangstaObjects}. Usually this is most useful if
#' \code{gangstaObjects} is set equal to results returned from
#' \code{getGangstas} or \code{subsetGangstas}.
#'
#' @param gangstaObjects List of \code{gangsta} objects to be queried.
#' @param attributeName Name of the attribute used for the query.
#' @param attributeValue All \code{gangsta} objects where the attribute
#'   (specified in \code{attributeName}) is equal to \code{attributeValue} will
#'   be returned.
#' @param gangstaNames Character vector containing the names (value stored in
#'   the \code{name} attribute) of the \code{gangsta} objects to be returned.
#' @param checkClass When set to a value other than "", checks to be sure all of
#'   the objects referenced in \code{gangstaNames} are of the class
#'   \code{checkClass}. Throws an error if not.
#' @param attribName Name of an attribute associated with 1 or more
#'   \code{gangsta} classes.
#' @return \code{gangstasExist} returns TRUE if all names in \code{gangstaNames}
#'   have an associated \code{gangsta} object in \code{gangstaObjects}.  If not,
#'   execution is halted and an error reported.  \code{getGangstaAttribute}
#'   returns a vector of attribute values from the objects in
#'   \code{gangstaObjects}.  All other functions return a \code{list} of
#'   \code{gangsta} objects.

#' @export
subsetGangstas = function(gangstaObjects, attributeName, attributeValue) {
  if(length(gangstaObjects) == 0) {
    returnVal = list()
  } else {
    if(attributeName == "class") {
      itMatches = sapply(gangstaObjects, is, class2 = attributeValue)
    } else {
      if(is.logical(attributeValue)) {
        itMatches = sapply(gangstaObjects, "[", name=attributeName)
        ## where itMatches is TRUE, return TRUE; where NULL or FALSE, return FALSE
        itMatches = sapply(itMatches, function(x) !is.null(x) && (!xor(x, attributeValue)))
      } else {
        itMatches = (sapply(gangstaObjects, "[", name=attributeName) == attributeValue)
      }
    }
    returnVal = gangstaObjects[itMatches]
  }
  return(returnVal)
}

#' @rdname subsetGangstas
#' @export
getGangstas = function(gangstaObjects, gangstaNames) {
  hits = lapply(gangstaNames, subsetGangstas, gangstaObjects = gangstaObjects, attributeName = gangstaAttributeName("name"))
  notFound = (sapply(hits, length) == 0)
  if(any(notFound)){
    stop(paste("Gangstas with the following names were requested but not found in the list of gangsta objects: ", paste0(gangstaNames[notFound], collapse = ", ")))
  }
  duplicates = (sapply(hits, length) > 1)
  if(any(duplicates)) {
    stop(paste0("More than one gangsta object exists with the following names ", paste0(gangstaNames[duplicates], collapse = ", ")))
  }
  return(unlist(hits, recursive = F))
}

#' @rdname subsetGangstas
#' @export
getGangstaAttribute = function(gangstaObjects, attribName) {
  sapply(gangstaObjects, "[[", attribName)
}

#' @rdname subsetGangstas
#' @export
gangstasExist = function(gangstaObjects, gangstaNames, checkClass = "") {
  ## The next line calls "getGangstas()", which will halt execution
  ## and throw and error if a gangsta doesn't exist.
  matchingGangstas = getGangstas(gangstaObjects, gangstaNames)
  if(any(checkClass!="")) {
    badClass = !sapply(matchingGangstas, is, checkClass)
    if(any(badClass)) {
      stop(paste0("Gangsta objects of type '", checkClass, "' are required for, but the following requested gangsta objects are not of that class: ", paste0(gangstaNames[badClass], collapse = ", ")))
    }
  }
  return(TRUE)
}

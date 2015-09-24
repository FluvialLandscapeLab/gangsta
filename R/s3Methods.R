print.gangsta = function(x) {
  classes = attr(x,"class")
  cat(classes[1], ":", x$name, "\n")
  cat("classes :", paste0(classes, collapse = ", "), "\n")
  otherVals = x[names(x)[names(x)!="name"]]
  mapply(function(nm, val) cat(nm, ":", paste0(addNamesToValues(val), collapse = "; "), "\n"), names(otherVals),otherVals)
}

addNamesToValues = function(x) {
  if(is.null(names(x))) {
    return(x)
  } else {
    return(paste(names(x), "=", x))
  }
}


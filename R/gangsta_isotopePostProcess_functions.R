exprsnPaste = function(...) {
  if(any(unlist(lapply(list(...), function(x) length(x) == 0)))) {
    result = character(0)
  } else {
    result = paste(...)
  }
  return(result)
}

exprsnPaste0 = function(...) {
  return(exprsnPaste(..., sep = ""))
}

makeGenericVars = function(prefixes, varTag, separator = ".") {
  return(exprsnPaste0(prefixes, separator, gangstaVarName(varTag)))
}

gangstaVarName = function(varTag) {
  return(getOption("gangsta.vars")[varTag])
}

makePoolStartMolVars = function(poolNames) {
  return(makeGenericVars(poolNames, "startSuffixPool"))
}

makePoolEndMolVars = function(poolNames) {
  return(makeGenericVars(poolNames, "endSuffixPool"))
}

makeTransferMolTransVars = function(transferNames) {
  return(makeGenericVars(transferNames,"transSuffix"))
}

print.gangsta = function(x) {
  classes = attr(x,"class")
  cat(classes[1], ":", x$name, "\n")
  cat("classes :", paste0(classes, collapse = ", "), "\n")
  otherVals = x[names(x)[names(x)!="name"]]
  mapply(function(nm, val) cat(nm, ":", paste0(addNamesToValues(val), collapse = "; "), "\n"), names(otherVals),otherVals)
}

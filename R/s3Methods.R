print.gangsta = function(x) {
  classes = attr(x,"class")
  cat(classes[1], ":", x$name, "\n")
  cat("classes :", paste0(classes, collapse = ", "), "\n")
  otherVals = x[names(x)[names(x)!="name"]]
  mapply(function(nm, val) cat(nm, ":", val, "\n"), names(otherVals),otherVals)
}


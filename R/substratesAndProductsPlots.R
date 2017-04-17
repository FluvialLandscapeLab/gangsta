inputs = function() {
  compoundOrder = c("Met", "Aut", "Het", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "O2")
  inputList = unlist(leakInListInput, recursive = F)
  inputList = inputList[sapply(inputList, length) != 0]
  nInputs = length(inputList)
  inputs = rep(NA, nInputs)
  for(i in 1:nInputs){
    inputs[i] = inputList[[i]]$additionalMols
    names(inputs)[i] = inputList[[i]]$compoundName
  }
  inputs = tapply(inputs, names(inputs), sum)
  inputNames = compoundOrder[compoundOrder %in% names(inputs)]
  inputs = inputs[match(inputNames, names(inputs))]
  return(inputs)
}

output = function(resultsList, gangstas){
  compoundOrder = c("Met", "Aut", "Het", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "O2")
  inputs = inputs()
  nIter = length(resultsList)
  compoundMols = resultsList[[nIter]]$compoundVals$final
  names(compoundMols) = row.names(resultsList[[nIter]]$compoundVals)
  compoundsToAdd = inputs[!names(inputs) %in% names(subsetGangstas(gangstas, "class", "compound"))]
  output = c(compoundMols, compoundsToAdd)
  output = output[output > 0]
  outputNames = compoundOrder[compoundOrder %in% names(output)]
  output = output[match(outputNames, names(output))]
  return(output)
}

substrProdPlot = function(resultsList, gangstas){
  input = inputs()
  output = output(resultsList, gangstas)
  inputDF = data.frame(type = "input", compound = names(input), umols = input)
  outputDF = data.frame(type = "output", compound = names(output), umols = output)
  inOutDF = rbind(inputDF, outputDF)

  backgroundCol = "white"
  textCol = "black"
  axisFontSize = 1.5

  substrAndProdPlot =
    ggplot2::ggplot(inOutDF, ggplot2::aes(x = type, y = umols, fill = compound)) +
    ggplot2::geom_bar(stat = "identity") +


    # ggplot2::expand_limits(y=c(0,0.1)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = backgroundCol, colour = backgroundCol),
                   panel.background = ggplot2::element_blank()) +
    # ggplot2::scale_x_discrete(
    #   inOutDF$type,
    #   labels = c("Substrates", "Products"),
    #   name = ggplot2::element_blank()
    # ) +
    # ggplot2::scale_y_continuous(
      # breaks = seq(0, 0.09, 0.03),
      # labels = as.character(seq(0, 0.09, 0.03)),
      # name = ""
    # ) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = textCol, size = ggplot2::rel(axisFontSize))) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = backgroundCol),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), colour = textCol  ),
      legend.position = "bottom") +
    ggplot2::scale_fill_manual(
      values = palette(rainbow(14)),
      labels = c("Met", "Aut", "Het", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "O2")
    )
    # ggplot2::theme(plot.margin = grid::unit(c(-0.60, 1, 0.02, 0.5), "lines"))
}

combineSubstrProdPlots = function() {
  resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"]
  resultNames = sort(resultNames)

  gangstaNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,8) == "gangstas"]
  gangstaNames = sort(gangstaNames)}

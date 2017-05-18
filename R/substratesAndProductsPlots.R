inputs = function(removeInfiniteCompounds = T, infCompoundNames = c("Ox", "Hx")) {
  compoundOrder = c("Met", "Aut", "Het", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "O2")
  inputList = unlist(leakInListInput, recursive = F)
  inputList = inputList[sapply(inputList, length) != 0]
  nInputs = length(inputList)
  inputs = rep(NA, nInputs)
  for(i in 1:nInputs){
    inputs[i] = inputList[[i]]$additionalMols
    names(inputs)[i] = inputList[[i]]$compoundName
  }
  # sum inputs for all time steps (e.g., if something is added at multiple time steps, sum it)
  inputs = tapply(inputs, names(inputs), sum)

  # order input compound names according to compoundOrder
  inputNames = compoundOrder[compoundOrder %in% names(inputs)]
  inputs = inputs[match(inputNames, names(inputs))]

  # create a vector of 0 for compound mols
  compounds = structure(rep(0, length(compoundOrder)), .Names = compoundOrder)

  # add mols for compounds which have input mols > 0
  compounds[match(names(inputs),names(compounds))] = inputs

  # remove infinite compounds
  if(removeInfiniteCompounds == T) {
    compounds = compounds[!names(compounds) %in% infCompoundNames]
  }
  return(compounds)
}

filteredInputs = function(resultsList, gangstas, removeInfiniteCompounds = T, infCompoundNames = c("Ox", "Hx")){
  #set order of compounds
  compoundOrder = c("Met", "Aut", "Het", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "O2")

  #get number of time steps
  nIter = length(resultsList)

  #this concise little nugget of code gets everything that was added during each
  #time step and sums it for the enitre model run.  Yes, Goff, you've created a
  #monster.
  inputs = apply(do.call(rbind, lapply(resultsList, function(ts) ts$leakInCompoundVals)), 2, sum)

  # order input compound names according to compoundOrder
  inputNames = compoundOrder[compoundOrder %in% names(inputs)]
  inputs = inputs[match(inputNames, names(inputs))]

  # create a vector of 0 for compound mols
  compounds = structure(rep(0, length(compoundOrder)), .Names = compoundOrder)

  # add mols for compounds which have input mols > 0
  compounds[match(names(inputs),names(compounds))] = inputs

  # remove infinite compounds
  if(removeInfiniteCompounds == T) {
    compounds = compounds[!names(compounds) %in% infCompoundNames]
  }
  return(compounds)
}

outputs = function(resultsList, gangstas, removeInfiniteCompounds = T, infCompoundNames = c("Ox", "Hx")){
  #set order of compounds
  compoundOrder = c("Met", "Aut", "Het", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2", "O2")

  #get number of time steps
  nIter = length(resultsList)

  # get outputs from results list and assign their names
  compoundMols = resultsList[[nIter]]$compoundVals$final
  names(compoundMols) = row.names(resultsList[[nIter]]$compoundVals)
  outputs = compoundMols

  # order output compound names according to compoundOrder
  outputNames = compoundOrder[compoundOrder %in% names(outputs)]
  outputs = outputs[match(outputNames, names(outputs))]

  # create a vector of 0 for compound mols
  compounds = structure(rep(0, length(compoundOrder)), .Names = compoundOrder)

  # add mols for compounds which have output mols > 0
  compounds[match(names(outputs),names(compounds))] = outputs

  # remove infinite compounds
  if(removeInfiniteCompounds == T) {
    compounds = compounds[!names(compounds) %in% infCompoundNames]
  }
  return(compounds)
}

substrProdPlot = function(resultsList, gangstas, withLegend = F){

  # decide if you want a legend
  if(withLegend ==F) {
    legendLocation = "none"
  }else{
    legendLocation = "top"
  }

  # get inputs and outputs and put them in a data frame
  input = inputs()
  filtInput = filteredInputs(resultsList, gangstas)
  output = outputs(resultsList, gangstas)
  inputDF = data.frame(type = "I", compound = names(input), umols = input)
  filtInputDF = data.frame(type = "R", compound = names(filtInput), umols = filtInput)
  outputDF = data.frame(type = "P", compound = names(output), umols = output)
  inOutDF = rbind(inputDF, filtInputDF, outputDF)

  # remove organisms from the data frame
  inOutDF = inOutDF[!(inOutDF$compound == "Met"|inOutDF$compound == "Aut"| inOutDF$compound == "Het"),]

  # this next bit of code is a major hack that I don't know how to get around.
  # We need a box for O2 in the reactants for the scenario wherein O2 is a
  # source/sink.  In the particular scenario that we want to put in the paper,
  # it is a CN model with O2, Ox, and Hx as source/sinks.  Well, O2 isn't listed
  # in the gangstas list as an infinite compound because the model is a carbon
  # and nitrogen model only.  Run the following line of code if you don't
  # believe me: getGangstaAttribute(subsetGangstas(gangstasCN_Ox.O2.Hx ,
  # "class", "compound"), "InfiniteCompound").  Okay, so I am setting the O2 box
  # to be of arbitrary size in the Reactants of the gangstasCN_Ox.O2.Hx model.
  if(identical(resultsList, resultsCN_Ox.O2.Hx)){
    inOutDF$umols[inOutDF$type == "R" & inOutDF$compound == "O2"] = 0.75
  }

  # make compound a factor with levels in the order that I want them in the plot
  # legend
  inOutDF$compound = factor(inOutDF$compound, levels = c("O2", "Met", "Aut", "Het", "DOM", "Ox", "CO2", "CH4", "SO4", "HS", "Hx", "NH4", "NO3", "N2"))

  colorVect = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

  # set background, text color, and font sizes of plot
  backgroundCol = "white"
  textCol = "black"
  axisFontSize = 1.5



  # make the plot
  substrAndProdPlot =
    ggplot2::ggplot(inOutDF, ggplot2::aes(x = type, y = umols, fill = compound)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::expand_limits(y=c(0,1.25)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = backgroundCol, colour = backgroundCol),
                   panel.background = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1.0, 0.5),
      labels = c("0", "0.5", "1.0"),
      name = ""
    ) +
    ggplot2::scale_fill_manual(
      values = colorVect
    )+
    ggplot2::scale_x_discrete(name="")+
    ggplot2::theme(axis.text = ggplot2::element_text(colour = textCol, size = ggplot2::rel(axisFontSize))) +
    # ggplot2::theme(
    #   axis.text.x = ggplot2::element_text(size = 20),
    #   axis.text.y = ggplot2::element_text(size = 20)
    # ) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = backgroundCol),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), colour = textCol  ),
      legend.position = legendLocation)
  return(substrAndProdPlot)
}

combineSubstrProdPlots = function(
  withLegend = F,
  fileIdx,
  axisFontSize = 1.5) {

  filePrefix = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\Figures\\SubAndProdPlots\\"
  fileName = makeFileName(fileID = fileIdx, filePrefix = filePrefix)

  resultNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,7) == "results"]
  resultNames = sort(resultNames)

  gangstaNames = ls(envir = .GlobalEnv)[substring(ls(envir =.GlobalEnv), 1,8) == "gangstas"]
  gangstaNames = sort(gangstaNames)

  plotList =
    base::lapply(
      1:length(resultNames),
      function(rNum)
        substrProdPlot(
          results =  get(resultNames[rNum], envir = .GlobalEnv),
          gangstas = get(gangstaNames[rNum], envir = .GlobalEnv),
          withLegend = withLegend
        )
    )
  ggplot2::ggsave(
    filename = fileName,
    plot =
      gridExtra::grid.arrange(
        grobs = plotList,
        ncol = 1, nrow = length(resultNames),
        heights = rep(2, length(resultNames))
      ),
    height = 8.5,
    width = 2,
    dpi = 600
  )
}





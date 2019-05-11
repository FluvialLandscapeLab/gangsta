isotopesToTrack = list(C = c(12,13))

RstC <- 0.0112372

initialIsotopicRatios = list(DIC_C = c("12" = 0.9,"13" = 0.1),
                             CH4_C = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)),
                             DOM_C = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)),
                             Acetate_C = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)),
                             Aut_C = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)),
                             Met_C = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)),
                             Het_C = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)),
                             Acetoclast_C = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))
isotopeLeakInList =  list(
  list(list()),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1)))),
  list(list(poolName = "Met_C", isotopicRatios = c("12" = 1-RstC/(RstC+1),"13" = RstC/(RstC+1))))
)

results <- isotopePostProcess(results = results,
                              gangstas = myGangstas,
                              isotopesToTrack = isotopesToTrack,
                              initialIsotopicRatios = initialIsotopicRatios,
                              isotopeLeakInList = isotopeLeakInList)

CO2_C13 = calculateDelVals(results = results,
                           poolName = "DIC_C",
                           AtomicWeight = "13",
                           RStd = RstC)
plotIsotopicComposition(results = results,
                        poolName = "DIC_C",
                        AtomicWeight = "13",
                        RStd = RstC)
plotAllIsotopicCompositions(results = results,
                            gangstas = myGangstas,
                            elementName = "C",
                            AtomicWeight = "13",
                            RStd = RstC)
plotTransfersIn(results, "DIC_C")
plotTransfersOut(results, "DIC_C")
plotTransfersIn(results, "DOM_C")
plotTransfersOut(results, "DOM_C")
plotTransfersIn(results, "CH4_C")
plotTransfersOut(results, "CH4_C")
plotAllTransfers(results, "DIC_C")
plotAllTransfers(results, "DOM_C")
plotAllTransfers(results, "Het_C")
plotAllTransfers(results, "Met_C")




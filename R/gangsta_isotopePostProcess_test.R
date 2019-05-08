isotopesToTrack = list(C = c(12,13))

RstC <- 0.0112372

initialIsotopicRatios = list(DIC_C = c("12" = 0.9,"13" = 0.1),
                             CH4_C = c("12" = 1-RstC,"13" = RstC),
                             DOM_C = c("12" = 1-RstC,"13" = RstC),
                             Acetate_C = c("12" = 1-RstC,"13" = RstC),
                             Aut_C = c("12" = 1-RstC,"13" = RstC),
                             Met_C = c("12" = 1-RstC,"13" = RstC),
                             Het_C = c("12" = 1-RstC,"13" = RstC),
                             Acetoclast_C = c("12" = 1-RstC,"13" = RstC))

results <- isotopePostProcess(results = results,
                              gangstas = myGangstas,
                              isotopesToTrack = isotopesToTrack,
                              initialIsotopicRatios = initialIsotopicRatios)

CO2_C13 = calculateDelVals(results = results,
                           poolName = "CO2_C",
                           AtomicWeight = "13",
                           RStd = RstC)
plotIsotopicComposition(results = results,
                        poolName = "CO2_C",
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



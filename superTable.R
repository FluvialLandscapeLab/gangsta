resultsList = gangstasCONSH_Ox.Hx

transfers.dat <-
  do.call(
    rbind.data.frame,
    subsetGangstas(resultsList, "class", "transformation")
    )

process.energies.dat <-
  as.data.frame(
    sapply(
      subsetGangstas(resultsList, "class", "process"),
      "[[",
      "energyTerm"
    )
  )
names(process.energies.dat) = "Affinity" #formerly called "energyTerm"
process.energies.dat$processName = row.names(process.energies.dat)

transfers.dat$organism <- substr(transfers.dat$processName, 1, 3)
transfers.dat$baseProcessName <- substr(transfers.dat$processName, 4, nchar(as.character(transfers.dat$processName)) )
transfers.dat$element <- substr(transfers.dat$name, nchar(as.character(transfers.dat$name)), nchar(as.character(transfers.dat$name)) )


transfers.dat <-
  merge(
    x = transfers.dat,
    y = process.energies.dat,
    by.x = "processName",
    by.y = "processName",
    all.x = T,
    all.y = F
  )

for(org in c("Het", "Aut", "Met")) {
  transfers.dat$from <- gsub(org, "Bio", transfers.dat$from)
  transfers.dat$to <- gsub(org, "Bio", transfers.dat$to)
}

transfers.dat$transferPools <- paste(transfers.dat$from, transfers.dat$to, sep = " --> ")

transfers.dat$multiToPools = F
for(i in 2:nrow(transfers.dat)) {
  if(
    i %in% grep("Assim", transfers.dat$processName) &
    transfers.dat$processName[i-1] ==  transfers.dat$processName[i] &
    transfers.dat$from[i-1] ==  transfers.dat$from[i] &
    transfers.dat$molarTerm[i-1] ==  transfers.dat$molarTerm[i] &
    transfers.dat$joulesToMolsRatio[i-1] ==  transfers.dat$joulesToMolsRatio[i]) {
    transfers.dat$multiToPools[i-1] = T
    transfers.dat$multiToPools[i] = T
  }
}

superTable <-
  plyr::ddply(
    transfers.dat,
    c("baseProcessName", "Affinity", "transferPools", "molarTerm", "joulesToMolsRatio", "multiToPools"),
    plyr::summarise,
    organisms = list(as.character(organism))
  )

superTableOutput <- data.frame(lapply(superTable, as.character), stringsAsFactors=FALSE)

write.csv(
  superTableOutput,
  file = "C:\\Users\\AnnMarie\\Dropbox\\GangstaShare\\gangstaManuscript\\ProcessesTable\\superTable.csv"
  )

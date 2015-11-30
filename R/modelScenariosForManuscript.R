runIt = function() {
  pdf("modelResults_20151130.pdf",
      onefile = T,
      paper = "USr",
      width = 10.5,
      height = 8
      )

  # Anaerobic carbon
  CNOS_Any("C")
  # Anaerobic N
  CNOS_Any("N")
  # Anaerobic C,N
  CNOS_Any(c("C", "N"))

  # Aerobic carbon
  CNOS_Any("C", "O2")
  # Aerobic N
  CNOS_Any("N", "O2")
  # Aerobic C,N
  CNOS_Any(c("C", "N"), "O2")

  # C,O
  CNOS_Any(c("C", "O"))
  # N,O
  CNOS_Any(c("N","O"))
  # C,N,O
  CNOS_Any(c("C", "N", "O"))
  # C,N,O,S
  CNOS_Any(c("C", "N", "O", "S"))

  dev.off()
}

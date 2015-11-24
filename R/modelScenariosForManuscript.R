runIt = function() {
  pdf("modelResults_20151124.pdf",
      onefile = T,
      paper = "USr",
      width = 10.5,
      height = 8
      )

  # Anaerobic carbon
  CNOS_Any("C", "Ox")
  # Anaerobic N
  CNOS_Any(c("N"), "Ox")
  # Anaerobic C,N
  CNOS_Any(c("C", "N"), "Ox")

  # Aerobic carbon
  CNOS_Any("C", c("O2", "Ox"))
  # Aerobic N
  CNOS_Any(c("N"), c("O2", "Ox"))
  # Aerobic C,N
  CNOS_Any(c("C", "N"), c("O2", "Ox"))

  # C,O
  CNOS_Any(c("C", "O"), c("Ox"))
  # C,N,O
  CNOS_Any(c("C", "N", "O"), "Ox")
  # C,N,O,S
  CNOS_Any(c("C", "N", "O", "S"), "Ox")

  dev.off()
}

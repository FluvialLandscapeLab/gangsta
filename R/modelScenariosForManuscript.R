runIt = function() {
  pdf("modelResults.pdf", onefile = T, paper = "USr", width = 10.5, height = 8)
  # Anaerobic carbon
  CNOS_Any("C", "Ox")
  # Anaerobic C,N
  CNOS_Any(c("C", "N"), "Ox")

  # Aerobic carbon
  CNOS_Any("C", c("O2", "Ox"))
  # Aerobic C,N
  CNOS_Any(c("C", "N"), c("O2", "Ox"))

  # C,N,O
  CNOS_Any(c("C", "N", "O"), "Ox")
  # C,N,O,S
  CNOS_Any(c("C", "N", "O", "S"), "Ox")

  dev.off()
}

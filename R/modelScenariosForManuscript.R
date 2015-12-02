runIt = function() {
  pdf("modelResults_20151130.pdf",
      onefile = T,
      paper = "USr",
      width = 10.5,
      height = 8
      )

  # Anaerobic carbon
  CNOSH_Any("C")
  # Anaerobic N
  CNOSH_Any("N")
  # Anaerobic C,N
  CNOSH_Any(c("C", "N"))

  # Aerobic carbon
  CNOSH_Any("C", "O2")
  # Aerobic N
  CNOSH_Any("N", "O2")
  # Aerobic C,N
  CNOSH_Any(c("C", "N"), "O2")

  # C,O
  CNOSH_Any(c("C", "O"))
  # N,O
  CNOSH_Any(c("N", "O"))
  # C,N,O
  CNOSH_Any(c("C", "N", "O"))
  # C,N,O,S
  CNOSH_Any(c("C", "N", "O", "S"))

  CNOSH_Any(c("O", "H"))

  dev.off()
}

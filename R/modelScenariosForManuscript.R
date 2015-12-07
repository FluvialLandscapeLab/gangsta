runIt = function() {
  pdf("modelResults_20151203.pdf",
      onefile = T,
      paper = "letter",
      width = 8,
      height = 10.5
      )

  # Anaerobic carbon
  CNOSH_Any("C")
  # Anaerobic C,N,S
  CNOSH_Any(c("C", "N", "S"))

  # Aerobic C,N,S
  CNOSH_Any(c("C", "N" ,"S"), "O2")

  CNOSH_Any(c("C", "N", "S", "O"))

  CNOSH_Any(c("O", "H"))

  dev.off()
}

computeAsymmetry <- function(asymmetry, wind, x, y, vx, vy, vh, r, rmw, lat) {
  # Circular symmetrical wind
  dir <- -(atan2(y, x) - pi / 2)

  if (lat >= 0) {
    dir <- dir - pi / 2
  } else {
    dir <- dir + pi / 2
  }

  dir[dir < 0] <- dir[dir < 0] + 2 * pi
  dir[dir > 2 * pi] <- dir[dir > 360] - 2 * pi

  windX <- wind * cos(dir)
  windY <- wind * sin(dir)

  # Moving wind
  stormDir <- -(atan2(vy, vx) - pi / 2)
  if (is.na(stormDir)) browser()
  if (stormDir < 0) {
    stormDir <- stormDir + 2 * pi
  }

  mWindX <- vh * cos(stormDir)
  mWindY <- vh * sin(stormDir)

  # Formula for asymmetry
  if (asymmetry == "Chen") {
    formula <- 3 * rmw**(3 / 2) * r**(3 / 2) / (rmw**3 + r**3 + rmw**(3 / 2) * r**(3 / 2))
  } else if (asymmetry == "Miyazaki") {
    formula <- exp(-r / 500 * pi)
  }

  # New total wind speed
  tWindX <- windX + formula * mWindX
  tWindY <- windY + formula * mWindY


  wind <- sqrt(tWindX**2 + tWindY**2)

  # New wind direction
  direction <- atan2(tWindY, tWindX) * 180 / pi
  direction[direction < 0] <- direction[direction < 0] + 360

  return(list(wind = round(wind, 3), direction = round(direction, 3)))
}


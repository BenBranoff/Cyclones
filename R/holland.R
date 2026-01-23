holland <- function(L, cens,tmpRas) {
  cent <- cens |> slice(L)
  tmpRas <- crop(tmpRas,ext(unlist(st_drop_geometry(cent[,c('xmin','xmax','ymin','ymax')]))))
  msw = round(cent$maxWind/1.94384)
  rmw <- round(cent$rmw/1000)
  maxpress= cent$maxpress/0.01
  minpress=cent$minpress/0.01
  lat=cent$centerY_shifted_geo
  rho <- 1.15 # air densiy
  f <- 2 * 7.29 * 10**(-5) * sin(lat) # Coriolis parameter
  b <- rho * exp(1) * msw**2 / (maxpress - minpress)

  crds <- terra::crds(tmpRas, na.rm = FALSE)
  x <- crds[, 1] - cent$centerX_shifted_geo
  y <- crds[, 2] - cent$centerY_shifted_geo
  # Computing distances to the eye of the storm in m
  distEye <- terra::distance(
    x = crds,
    y = st_coordinates(cent|>st_transform(4326)|>st_shift_longitude()),
    lonlat = TRUE
  )* 0.001


  vr <- distEye
  vr <- sqrt(b / rho * (rmw / distEye)**b * (maxpress - minpress) * exp(-(rmw / distEye)**b) + (distEye * f / 2)**2) - distEye* f / 2
  tmpRasA <- computeAsymmetry("Chen",vr,x,y,cent$vxDeg_geo, cent$vyDeg_geo,
                              cent$stormSpeed_geo,
                              distEye, rmw, lat)
  #if (L>=79) browser()
  dist <- sqrt(x * x + y * y)
  tmpRasA$wind[dist > 2.5] <- NA
  tmpRasA$wind <- round(tmpRasA$wind,3)
  terra::values(tmpRas) <- tmpRasA$wind
  terra::time(tmpRas) <- cent$date
  #comps <- compare_winds(rasts=tmpRas,shape=data.frame(name=cent$name,date=cent$date))
  #comps$rsource="Holland"
  cat(paste0("\rCalculating Holland wind field: %",round(100*L/nrow(cens),1)))
  return(tmpRas)
}

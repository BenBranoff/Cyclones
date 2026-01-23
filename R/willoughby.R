#' @importFrom terra time
willoughby <- function(L,cens,tmpRas) {
  cent <- cens |> slice(L)
  tmpRas <- crop(tmpRas,ext(unlist(st_drop_geometry(cent[,c('xmin','xmax','ymin','ymax')]))))
  msw = round(cent$maxWind/1.94384)
  rmw <- round(cent$rmw/1000)
  lat=cent$centerY_shifted_geo
  crds <- terra::crds(tmpRas, na.rm = FALSE)
  x <- crds[, 1] - cent$centerX_shifted_geo
  y <- crds[, 2] - cent$centerY_shifted_geo
  # Computing distances to the eye of the storm in m
  distEye <- terra::distance(
    x = crds,
    y = st_coordinates(cent|>st_transform(4326)|>st_shift_longitude()),
    lonlat = TRUE
  )* 0.001

  x1 <- 287.6 - 1.942 * msw + 7.799 * log(rmw) + 1.819 * abs(lat)
  x2 <- 25
  a <- 0.5913 + 0.0029 * msw - 0.1361 * log(rmw) - 0.0042 * abs(lat)
  n <- 2.1340 + 0.0077 * msw - 0.4522 * log(rmw) - 0.0038 * abs(lat)

  vr <- distEye
  vr[distEye >= rmw] <- msw * ((1 - a) * exp(-abs((distEye[distEye >= rmw] - rmw) / x1)) + a * exp(-abs(distEye[distEye >= rmw] - rmw) / x2))
  vr[distEye < rmw] <- msw * abs((distEye[distEye < rmw] / rmw)^n)
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
  #comps$rsource="Willoughby"
  cat(paste0("\rCalculating Willoughby wind field: %",round(100*L/nrow(cens),1)))
  return(tmpRas)
}

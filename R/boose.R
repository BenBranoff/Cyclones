#' @importFrom terra ext
#' @importFrom dplyr nth
boose <- function(L,cens,tmpRas) {
  vx <- cens |> pull(vxDeg_geo)|>nth(L)
  vy <- cens |> pull(vyDeg_geo)|>nth(L)
  vh <- cens |> pull(stormSpeed_geo)|>nth(L)
  cent <- cens |> slice(L)
  tmpRas <- crop(tmpRas,ext(unlist(st_drop_geometry(cent[,c('xmin','xmax','ymin','ymax')]))))
  #custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",
  #                  st_coordinates(line|>st_transform(4326))[,"X"],
  #                  " +lat_0=",
  #                  st_coordinates(line|>st_transform(4326))[,"Y"])
  rho <- 1 # air densiy
  msw = cent$maxWind/1.94384
  maxpress= cent$maxpress/0.01
  minpress=cent$minpress/0.01
  rmw <- cent$rmw/1000
  land <- rnaturalearth::ne_countries() |> st_cast("POLYGON")|>st_transform(crs(tmpRas))|>rasterize(y=tmpRas,field=1,background=0)
  # Computing coordinates of raster
  crds <- terra::crds(tmpRas, na.rm = FALSE)
  x <- crds[, 1] - cent$centerX_shifted_geo
  y <- crds[, 2] - cent$centerY_shifted_geo
  # Computing distances to the eye of the storm in m
  distEye <- terra::distance(
    x = crds,
    y = st_coordinates(cent|>st_transform(4326)|>st_shift_longitude()),
    lonlat = TRUE
  )* 0.001
  b <- rho * exp(1) * msw**2 / (maxpress-minpress)
  vr <- distEye
  vr <- sqrt((rmw / vr)**b * exp(1 - (rmw / vr)**b))
  if (cent$centerY_shifted_geo >= 0) {
    # Northern Hemisphere, t is clockwise
    angle <- atan2(vy, vx) - atan2(y, x)
  } else {
    # Southern Hemisphere, t is counterclockwise
    angle <- atan2(y, x) - atan2(vy, vx)
  }
  landIntersect <- extract(land,crds)
  vr[landIntersect == 1] <- 0.8 * (msw - (1 - sin(angle[landIntersect == 1])) * vh[landIntersect == 1] / 2) * vr[landIntersect == 1]
  vr[landIntersect == 0] <- (msw - (1 - sin(angle[landIntersect == 0])) * vh / 2) * vr[landIntersect == 0]
  vr <- round(vr,3)
  dist <- sqrt(x * x + y * y)
  vr[dist > 2.5] <- NA
  terra::values(tmpRas) <- vr
  comps <- compare_winds(rasts=tmpRas,shape=data.frame(name=cent$name,date=cent$date,meth="Boose"))
  if (!is.null(comps$rsamps)) comps$rsource="Boose"
  cat(paste0("\rCalculating Boose wind field: %",round(100*L/nrow(cens),1)))
  return(list(rasters=tmpRas,comps=comps))
}

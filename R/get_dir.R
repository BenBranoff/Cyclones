#' @importFrom terra vect project crs
#' @importFrom geosphere bearing
#' @keywords internal
get_dir <- function(r,C){
  r_points <- as.data.frame(r,xy=TRUE)
  r_vec <- vect(r_points, crs=crs(r), geom=c("x", "y"))
  r_vec <- project(r_vec, "epsg:4326")
  ###  from geosphere package bearing function
  r_vec$bearing <- bearing(r_vec,C |> st_transform(4326) |>st_coordinates())
  r_vec$WD <- (r_vec$bearing+ 90 - 180)%% 360
  r_vec <- project(r_vec,crs(r))
  WD <- rasterize(r_vec,r,field="WD")
  WD
}

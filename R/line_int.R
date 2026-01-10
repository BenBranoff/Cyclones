#' @importFrom terra rast ext crop rasterize interpolate mask
#' @importFrom fields Tps
line_int <- function(line,center,ex,s_res){
  ###  create empty raster, with the extent of the current wind field
  r <- rast(ext=ex, resolution=s_res,
            crs=crs(line))
  Beta <- 1.0036+0.0173*as.numeric(max(c(line$kts,line$maxWind))+
                                     0.0313*log(min(line$dist_m_mean,na.rm=T))+
                                     0.0087*st_coordinates(center%>%st_transform((4326)))[,"Y"])
  ##  now the pressure at each swath extent

  line <- line %>% ungroup()%>%
    mutate(P=minpress+((maxpress-minpress)*exp(-(min(dist_m_mean[kts>0])/dist_m_mean)^Beta)))

  ### crop the line to the outer extent of the storm
  outer=line %>% filter(quad=="ROCI")
  ###   Using the original ROCI can create extreme and unnatural shifts in wind Velocity with the thin spline method
  ###   to avoid this, bump the 0 velocity line out a bit further
  ###  we will still zero out the original ROCI later, this is only for the thine spline interpolation
  line <- bind_rows(line,
                    st_buffer(outer%>%st_cast("POLYGON"),outer$dist_m_mean*.5) %>%
                      mutate(quad="ROCI2")%>%st_cast("LINESTRING"))

  r <- crop(r,line%>%filter(quad=="ROCI2"))
  ####  remove the eye and the outer storm limits
  ###  we can set those to zero later
  ###  rasterize the remaining swaths, one for velocity and one for pressure
  rv <-  rasterize(line%>%filter(!line$quad %in%c("ROCI","eye")), r, "kts",touches=T,fun="max")
  rP <-  rasterize(line%>%filter(!line$quad %in%c("ROCI","eye")), r, "P",touches=T)
  ###  convert the xyz values to a dataframe for thin spline interpolation
  xyv <- as.data.frame(rv, xy=T,na.rm=F)
  xyP <- as.data.frame(rP, xy=T,na.rm=F)
  options(warn = -1)
  ###  fit the thin spline interpolation model
  #boolFalse<-F
  #atry=200
  #while(boolFalse==F)
  #{
  #  tryCatch({
  #tps_v <- fastTps(xyv[,1:2], aRange=1500,xyv[,3],lon.lat = TRUE)#fields::Tps(xyv[,1:2], xyv[,3])
  #tps_P <-  fastTps(xyP[,1:2], aRange=1500,xyP[,3],lon.lat = TRUE)#fields::Tps(xyP[,1:2], xyP[,3]);
  #    boolFalse<-T
  #  },error=function(e){atry=atry+100
  #  },finally={})
  #}
  tps_v <- Tps(xyv[,1:2], xyv[,3])
  tps_P <-  Tps(xyP[,1:2], xyP[,3])
  ###  use the model to predict the unknown wind speeds
  p_v <-interpolate(r, tps_v)
  ###  in case the model predicts higher than recorded wind speeds, set them to the maximum known
  p_v[p_v>max(line$kts,na.rm=T)] <- max(line$kts,na.rm=T)
  p_P <- interpolate(r, tps_P)
  ####  do the same for the minimum values
  p_P[p_P<min(line$minpress,na.rm=T)] <- min(line$minpress,na.rm=T)
  p_P[p_P>max(line$maxpress,na.rm=T)] <- max(line$maxpress,na.rm=T)
  p_v[p_v<0] <- 0
  ###  now set the outer limits of the storm to 0
  p_v <- mask(p_v,outer %>% st_cast("POLYGON"))
  ###  and if there is an eye, set it to 0 as well
  if ("eye" %in% unique(line$quad)){
    p_v <- mask(p_v,line %>% filter(quad=="eye")%>%st_cast("POLYGON"),inverse=T)
    names(p_v) <- "lyr.1"
  }
  windir <- get_dir(p_v,center)
  ##  to calculate power:
  ###  the density of air is a conversion function of the pressure
  ##  assuming the temperature is around 25 C or 298K
  dens <- p_P*100/(287.058*298)
  ##  power is then a function of the density and the pressure
  ##  https://www.e-education.psu.edu/emsc297/node/649
  ###  here, wind velocity in knots is converted to m/s
  ##  and everything is converted from Watts to kWh by multiplying by the time duration
  ## we are calculating the power on 1 square meter of a surface
  kW <- (0.5*dens*1*(p_v/1.94384)^3)/1000
  kW[is.na(kW)] <- 0
  list(power=kW,vel=p_v,dir=windir)
}

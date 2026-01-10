#' @importFrom terra rast ext extend approximate setValues shift resample app allNA ifel time
move_and_interp <- function(L,lines,centers,s_res=5000){
  if (L==length(unique(lines$date))){return(NULL)}
  lines$date <- as.POSIXct(lines$date,tz="GMT")
  ###  find the outer bounds of the storm
  outerext<- lines %>%filter(quad=="track") %>% summarise() %>%st_buffer(max(lines$dist_m_mean[linestrings$quad=="ROCI"],na.rm=T),nQuadSegs=2)
  lne <- lines %>% filter(quad!="track")
   #lne <<- lines %>%
    #  arrange(date,kts) %>%
    #  group_by(date) %>%
    #  mutate(l = seq(1,n())) %>%
    # ungroup() %>%
   ##  get the time difference to the next timestep
    #  mutate(dt=as.numeric(difftime(date,diag(sapply(l, function(n) dplyr::lag(date, n=n))))))
  ###  get the date stamps of the lines
  dts <- unique(lne$date)#[1:(length(unique(lne$date)))]
  d1 <- unique(dts)[L]
  d2 <- unique(dts)[which(unique(dts)==d1)+1]
  timestep <- as.numeric(difftime(d2,d1,units="hours"))
  if (!is.na(timestep)&timestep>.25){
    line1 <- lne %>%
      filter(date==d1)
    line2 <- lne %>%
      filter(date==d2)
    if(d1==unique(dts)[length(unique(dts))]){end=T}else{end=F}
  }else{
    stop("insufficient timestep in line segments")
  }
  ###  create custom crs centered on current track segment
  ##  this reduces geometry calculation errors from using a general global crs
  custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",
                    st_coordinates(centers[centers$ISO_TIME==d2,] %>%st_transform(4326))[,"X"],
                    " +lat_0=",
                    st_coordinates(centers[centers$ISO_TIME==d2,]%>%st_transform(4326))[,"Y"])
  ###  calculate the delta x and delta y between the two locations
  ###  we will use this to align the two dates on top of eachother so we can then interpolate between them
  coord_dif <- st_coordinates(centers[centers$ISO_TIME==d2,] %>% st_transform(custCRS))-
    st_coordinates(centers[centers$ISO_TIME==d1,]%>% st_transform(custCRS))
  ## get all of the dates and calculate the time difference between each
  ##  this is important for calculating energy dissipation, which is a function of the time
  dates <- centers$ISO_TIME[centers$ISO_TIME<=d2&centers$ISO_TIME>=d1]
  dt <- unique(difftime(dates, lag(dates),unit="hours"))
  dt <- as.numeric(dt[!is.na(dt)])
  ###  now shift the second line segment to be centered on the first
  ##  this will allow us to interpolate raster representations between the two times
  line2.2 <- st_as_sfc(line2 %>% st_transform(custCRS))-matrix(data=coord_dif,ncol=2)
  line2.2 <- st_sf(geom=line2.2)
  st_crs(line2.2) <- custCRS
  line2.2 <- line2.2 %>% st_transform(st_crs(line1))  %>%
    st_as_sf() %>%
    rename(geometry=geom)
  line2.2<- cbind(line2.2,st_drop_geometry(line2))
  ##  create the wind velocity and power fields for the two segments,
  ##  which are now centered over the same geographical space
  rs1 <- line_int(line1,centers[centers$ISO_TIME==d1,],ext(outerext),s_res)
  ###  include the time duration in the power calculation
  rs1$kW = rs1$kW*dt
  rs2 <- line_int(line2.2,centers[centers$ISO_TIME==d1,],ext(outerext),s_res)
  rs2$kW = rs2$kW*dt
  ####  get the comparison with the shapefile used to construct the windfield
  ####   this is for quality control
  comps <- compare_winds(rasts=list(rs1,rs2),shapes=list(line1,line2.2))
  ###  stretch the velocity and power fields so they have the same spatial extent
  ###  this is necessary to stack and interpolate between them
  rs1$vel <- extend(rs1$vel,rs2$vel)
  rs2$vel <- extend(rs2$vel,rs1$vel)
  rs1$power <- extend(rs1$power,rs2$power)
  rs2$power <- extend(rs2$power,rs1$power)
  rs1$dir <- extend(rs1$dir,rs2$dir)
  rs2$dir <- extend(rs2$dir,rs1$dir)
  ###   create a list of empty rasters to hold all of the wind fields
  ###   that represent the interval times between the two segments
  vlist <- list()
  plist <- list()
  dlist <- list()
  vlist[[1]] <- rs1$vel
  plist[[1]] <- rs1$power
  dlist[[1]] <- rs1$dir
  for (i in 2:(length(dates)-1)){
    vlist[[i]] <- setValues(rs1[[2]],NA)
    plist[[i]] <- setValues(rs1[[1]],NA)
    dlist[[i]] <- setValues(rs1[[3]],NA)
  }
  ##  now add in the end raster with known values from the second line segment
  vlist[[length(vlist)+1]] <- rs2[[2]]
  plist[[length(plist)+1]] <- rs2[[1]]
  dlist[[length(dlist)+1]] <- rs2[[3]]
  ###  load them to be processed as a raster
  vs <- rast(vlist)
  ps <- rast(plist)
  ds <- rast(dlist)
  ###  inerpolate (linearly) the missing values between the two end "points"
  ##  we assume here that velocity and power transitions are linear from one time point to the next
  ##  in this case, its a time difference of usually three hours
  vs <- approximate(vs)
  ps <- approximate(ps)
  ds <- approximate(ds)
  ###  these rasters now hold the interpolated velocity and power values for times between the two end points
  ##  but their geographical information is identical
  ##  they now need to be moved to their respective geographic locations
  ###  create a new empty list to hold the shifted rasters
  Vlist <- list()
  Plist <- list()
  Dlist <- list()
  ##  set the first raster to that of the calculated field from the first segment, which has the correct geographic representation
  ##  the second segment will not be included because it is included on the next iteration as the first segment
  ##  it was only used here as an endpoint for the approximation
  Vlist[[1]] <- vs[[1]]
  Plist[[1]] <- ps[[1]]
  Dlist[[1]] <- ds[[1]]
  ##  for the remaining rasters, shift them in space so they are centered on the points corresponding to their
  ##  respective timestamps
  for (i in 2:(length(dates)-1)){
    d1 <- dates[1]
    d2 <- dates[i]
    ###  calculate the dx and dy between the two center points
    coord_dif <- st_coordinates(centers[centers$ISO_TIME==d2,])-
      st_coordinates(centers[centers$ISO_TIME==d1,])
    ##  shift the rasters accordingly
    v <- shift(vs[[i]],coord_dif[1], coord_dif[2])
    p <- shift(ps[[i]],coord_dif[1], coord_dif[2])
    d <- shift(ds[[i]],coord_dif[1], coord_dif[2])
    ## add the shifted rasters to the list
    Vlist[[length(Vlist)+1]] <- v
    Plist[[length(Plist)+1]] <- p
    Dlist[[length(Dlist)+1]] <- d
  }
  ###  if this is the end segment, use the last raster from the first list
  if (end==T){
    Vlist[[length(Vlist)+1]] <- vs[[length(vs)]]
    Plist[[length(Plist)+1]] <- ps[[legth(ps)]]
    Dlist[[length(Dlist)+1]] <- ds[[legth(ds)]]
  }
  ###  set the geometry back to the original rasters, to ensure all are equal before further processing
  Plist <- lapply(Plist, function(x) resample(x,rs1$power))
  ###  set NA values to 0, unless all every instance of that pixel is NA, then it remains NA
  a <- allNA(rast(Plist))
  Plist <- ifel(is.na(rast(Plist)), ifel(a, NA, 0), rast(Plist))
  Vlist <- lapply(Vlist, function(x) resample(x,rs1$vel))
  a <- allNA(rast(Vlist))
  Vlist <- ifel(is.na(rast(Vlist)), ifel(a, NA, 0), rast(Vlist))
  Dlist <- lapply(Dlist, function(x) resample(x,rs1$dir))
  a <- allNA(rast(Dlist))
  Dlist <- ifel(is.na(rast(Dlist)), ifel(a, NA, 0), rast(Dlist))
  ###  for parallel processing, they need to be returned in lists
  V <- Vlist;P=Plist;D=Dlist
  V <- extend(V,rast(ext= outerext, resolution=s_res,crs=crs(lines)))
  if (end==T){
    terra::time(V) <- terra::time(P) <- terra::time(D) <- dates
  }else{
    terra::time(V) <- terra::time(P) <- terra::time(D) <- dates[1:(length(dates)-1)]
  }
  cat(paste("\rBuilding rasters: %",round(100*L/ length(unique(lines$date)),1)))
  Sys.sleep(0.01)
  list(Vel=V,Power=P,Dir=D,comps=data.frame(comps))
}

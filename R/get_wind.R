#' @importFrom sf st_coordinates st_as_sfc st_sf st_crs st_shift_longitude
#' @importFrom terra rast res ext extend approximate setValues shift resample app allNA ifel time
get_wind <- function(lines,s_res=20000,t_res=60,methods=NULL){
  if (is.null(dim(lines))){
    warning("More than one storm in data. Only first will be used. Use lapply or a parallel equivalent to repeat for multiple storms")
    lines <- lines[[1]]
  }
  winds <- list()
  centers = interp_track(lines,t_res)
  if(methods=="all") methods=c("TPS","Willoughby","Boose","Holland")

  if ("TPS" %in% methods||is.null(methods)){
    if (!"linestrings" %in% lines$extent_type) stop("wrong geometry type. TPS accepts only linestrings.")

    msw <- lapply(seq_along(unique(lines$date)),TPS_int,lines, centers, s_res)
    msw <- rast(msw)
    names(msw) <- rep(paste("MSW",unique(lines$name),unique(format(lines$date,"%Y")),"TPS",sep="_"),terra::nlyr(msw))
    msw <- wrap(msw)
    winds=append(winds,list(msw))
    names(winds)[length(winds)] <- paste("MSW",unique(lines$name),unique(format(lines$date,"%Y")),"TPS",sep="_")
    cat("\n")
  }
  if (any(grepl("Willoughby|Boose|Holland",methods))){
    rTemp <- rast(ext=lines |> ext(), resolution=s_res,
                  crs=crs(lines))
    rTemp <- project(rTemp,"epsg:4326")
    rTemp <- rast(ext=st_shift_longitude(st_transform(lines,4326))|> ext(), resolution=0.166667,#res(rTemp),
                  crs=crs(st_shift_longitude(st_transform(lines,4326))))
    centers <- centers |> left_join(lines |>filter(!quad %in% c("track","track_point"))|>
                                      group_by(date) |>
                                      summarise(name=unique(name),rmw = min(dist_m_min[dist_m_min>0],na.rm=TRUE),maxWind=max(maxWind,na.rm=TRUE),
                                                minpress=min(minpress,na.rm=TRUE),maxpress=max(maxpress,na.rm=TRUE),
                                                roci =max(dist_m_mean,na.rm=TRUE),
                                                geoext=list(st_bbox(st_shift_longitude(st_transform(geometry[quad=="ROCI"],crs(rTemp))))),
                                                xmin=unlist(geoext)[1],ymin=unlist(geoext)[2],xmax=unlist(geoext)[3],ymax=unlist(geoext)[4])|>
                                      select(-geoext)|>st_drop_geometry(),by=join_by(date)) |>
      ungroup() |>
      mutate(name=unique(name[!is.na(name)]),across(c(rmw:ymax),zoo::na.approx),
             dt=as.numeric(difftime(lead(date),date,units="hours")))|>
      ###  get the shifted coords necessary for the equation methods
      st_transform(4326) |>
      st_shift_longitude()
    centers$centerX_shifted_geo <- st_coordinates(centers)[,1]
    centers$centerY_shifted_geo <- st_coordinates(centers)[,2]
    centers$stormSpeed_geo <- terra::distance(
      x=vect(cbind(centers$centerX_shifted_geo, centers$centerY_shifted_geo),crs="epsg:4326"),
      y=vect(cbind(dplyr::lead(centers$centerX_shifted_geo), dplyr::lead(centers$centerY_shifted_geo)),crs="epsg:4326"),
      pairwise=TRUE) * (0.001 / centers$dt) / 3.6
    centers <- centers |>
      mutate(vxDeg_geo=(dplyr::lead(as.numeric(centerX_shifted_geo)) - as.numeric(centerX_shifted_geo)) / dt,
             vyDeg_geo=(dplyr::lead(as.numeric(centerY_shifted_geo)) - as.numeric(centerY_shifted_geo)) / dt) |>
      st_transform(crs(rTemp)) |> slice(-n())
    for (meth in methods[grep("Willoughby|Boose|Holland",methods)]){
      ###  to convert projected resolution to geo resolution, make a transformed projected template raster

      if (meth=="Boose"){
        msw <- lapply(seq_along(unique(centers$date)),boose,centers,rTemp)
      }else if (meth=="Willoughby"){
        msw <- lapply(seq_along(unique(centers$date)),willoughby,centers,rTemp)
      } else if (meth=="Holland"){
        msw <- lapply(seq_along(unique(centers$date)),holland,centers,rTemp)
      }
      #msw <- do.call(rbind,lapply(msw,'[[',1))
      #compare <- do.call(rbind,lapply(msw,'[[',2))
      #comps <- append(comps,compare)
      # Applying focal function to smooth results
      # from stormR
      nbgmod <- nls(y ~ SSasymp(x, Asym,R0, lrc),data=data.frame(y=c(59,11,5,3),x=c(1000,4500,9000,20000)))
      nbg <- predict(nbgmod,newdata=c(x=s_res))
      msw <- lapply(msw,function(x,y) {y=rast(ext=y |> ext(), resolution=s_res, crs=crs(y));project(x,y)},y=lines)
      msw <- rast(msw)
      names(msw) <- rep(paste("MSW",unique(lines$name),unique(format(lines$date,"%Y")),meth,sep="_"),terra::nlyr(msw))
      #msw <- max(rast(msw), na.rm = TRUE)
      msw <- terra::focal(msw, w = matrix(1, nbg, nbg), mean, na.rm = TRUE, pad = TRUE)
      msw <- wrap(msw)
      winds=append(winds,list(msw))
      names(winds)[length(winds)] <- paste("MSW",unique(lines$name),unique(format(lines$date,"%Y")),meth,sep="_")
      cat("\n")
    }
  }
  return(winds)
}



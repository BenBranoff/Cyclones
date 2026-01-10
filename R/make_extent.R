#' @importFrom dplyr last rowwise add_rownames join_by pull
#' @importFrom sf st_sfc st_point st_transform st_polygon st_linestring st_as_sf st_buffer st_cast st_set_geometry
#' @keywords internal
make_extent <- function(L,mods,tracks){
  tracks <- tracks |>
    mutate(ISO_TIME =as.POSIXct(ISO_TIME,format="%Y-%m-%d %H:%M:%S",tz="GMT"),
           MONTH  = format(ISO_TIME,"%m"))
  track <- tracks |> slice(L)
  if (nrow(track)>1){stop("more than one row in data. Input only one row at a time. Use lapply or a parallel equivalent to repeat for multiple rows")}
  options(dplyr.summarise.inform = FALSE)
  mod <- mods$modquads
  emod <- mods$emod
  pocimod <- mods$pocimod
  rocimod <- mods$rocimod
  minpresss <- mods$minpress
  id=track$SID
  name=track$NAME
  year=track$SEASON
  date <- track$ISO_TIME
  basin=track$BASIN
  ###  north american (NA) basins may return NA instead of  "NA"
  if(is.na(basin)){basin="NA"}
  X <- as.numeric(track$LON)
  Y <- as.numeric(track$LAT)
  ###  create a custom CRS centered on the track. This helps preserve true distances and removes error
  ##  associated with crs warping
  custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",X," +lat_0=",Y)
  ###  create another CRS for the whole set of tracks for this storm
  tracks$LAT <- as.numeric(tracks$LAT)
  tracks$LON <- as.numeric(tracks$LON)
  custCRSall <- paste0("+proj=lcc +lat_0=",round(mean(tracks$LAT,na.rm=TRUE))," +lon_0=",round(mean(tracks$LON,na.rm=TRUE))," +lat_1=",round((max(tracks$LAT,na.rm=TRUE)-min(tracks$LAT,na.rm=TRUE))/6),
                      " +lat_2=",round(5*(max(tracks$LAT,na.rm=TRUE)-min(tracks$LAT,na.rm=TRUE))/6)," +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
  ## get the center point of the current track
  center <- st_sfc(st_point(c(X,Y),dim="XY"),crs = 4326) |>
    st_transform(custCRS)
  ###  take relevant information from the IBTrACS data
  Cat <- as.numeric(track$USA_SSHS)
  maxwind <- as.numeric(track$USA_WIND)
  ###  in rare, mostly minor storms, the maxwind is missing
  ##   assign this based on the category
  if(is.na(maxwind)){
    maxwind=ifelse(Cat==5,137,ifelse(Cat==4,113,ifelse(Cat==3,96,ifelse(Cat==2,83,ifelse(Cat==1,64,ifelse(Cat==0,34,ifelse(Cat==-1,24,ifelse(Cat==-2,24,ifelse(Cat==-3,14,10)))))))))
  }
  maxwinddist <- as.numeric(track$USA_RMW)
  if (length(maxwinddist)==0){maxwinddist <- NA}
  minpress <- as.numeric(track$USA_PRES)
  if (length(minpress)==0){minpress <- NA}
  #### in only a very few cases, the minimum pressure is missing,
  ####  take the average for the storms of the same size class
  if (is.na(minpress)){
    minpress=minpresss$USA_PRES[minpresss$USA_SSHS==Cat&minpresss$MONTH==track$MONTH&minpresss$BASIN==basin]
    if (length(minpress)==0){
      minpress=mean(minpresss$USA_PRES[minpresss$MONTH==track$MONTH&minpresss$BASIN==basin],na.rm=TRUE)
    }
  }
  maxpress <- as.numeric(track$USA_POCI)
  if (length(maxpress)==0){maxpress=NA}
  ### there are more cases of  missing outer pressure, so model this
  if (is.na(maxpress)){
    if (length(pocimod$model[pocimod$USA_SSHS==Cat&pocimod$BASIN==basin])>0){
      maxpress <- minpress+predict(pocimod$model[pocimod$USA_SSHS==Cat&pocimod$BASIN==basin][[1]],newdata=data.frame(USA_WIND=maxwind))
    }else{
      maxpress <- minpress
    }
  }
  ROCId <- as.numeric(track$USA_ROCI)
  if (length( ROCId)==0){ ROCId=NA}
  EYE <- as.numeric(track$USA_EYE)
  if (length(EYE)==0){ EYE=NA}
  swathpoints <- c(NE34 = track$USA_R34_NE,SE34 = track$USA_R34_SE,NW34 = track$USA_R34_NW,SW34 = track$USA_R34_SW,
                   NE50 = track$USA_R50_NE,SE50 = track$USA_R50_SE,NW50 = track$USA_R50_NW,SW50 = track$USA_R50_SW,
                   NE64 = track$USA_R64_NE,SE64 = track$USA_R64_SE,NW64 = track$USA_R64_NW,SW64 = track$USA_R64_SW)
  swathpoints <- sapply(swathpoints,as.numeric)
  swathpoints[is.na(swathpoints)] <-0
  ###  store any native given information for later marking native versus modeled info
  swathmodeled <- swathpoints
  if (sum(swathpoints,na.rm=T)==0){
    ###  if there are no quadrant data, create a windfield based on the model from the minimum pressure
    ##  if the max wind and the min pressure are missing, there's not enough information
    ##  not much we can do
    swathpoints <- data.frame(speed=rep(c(34,50,64,maxwind),4),quad=rep(c("NE","SE","SW","NW"),each=4),USA_SSHS=Cat,BASIN=basin,source="modeled")
    swathpoints <- swathpoints[swathpoints$speed<=maxwind,]
    # models could only be fit for category 1 storms and above
    if (Cat>0){
      swathpoints <- swathpoints |>
        group_nest(quad,USA_SSHS,BASIN) |>
        left_join(mod|>select(quad,USA_SSHS,BASIN,model_asymp),by=c("USA_SSHS","quad","BASIN")) |>
        mutate(dist=map2(model_asymp,data,predict))|>
        select(-model_asymp) |>
        tidyr::unnest(c(data,dist))|>
        mutate(dist=round(dist))
      ##  for smaller storms, use the cat 1 storm model
    }else{
      swathpoints <- swathpoints |>
        group_nest(quad,BASIN) |>
        left_join(mod|>filter(USA_SSHS==1)|>select(quad,model_asymp,BASIN),by=c("quad","BASIN")) |>
        mutate(dist=map2(model_asymp,data,predict))|>
        select(-model_asymp) |>
        tidyr::unnest(c(data,dist))|>
        mutate(dist=round(dist))
    }
    swathpoints <- setNames(swathpoints$dist,paste0(swathpoints$quad,swathpoints$speed))
    swathpoints[swathpoints<10] <-10
  }else{
    ##  otherwise, use the given windfield, plus the maximum wind
    ##  if no RMW is given,we assume the maximum wind mirrors the extent of the maximum swath, but with smaller distance
    repl <- swathpoints[grep(max(as.numeric(gsub("NE|SE|SW|NW","",names(swathpoints[swathpoints>0]))),na.rm=T),names(swathpoints))]
    repl <- repl[repl>0]
    if (is.na(maxwinddist)){
      if (Cat>0){
        repl <- data.frame(speed=as.numeric(gsub("NE|SE|SW|NW","",names(repl))),
                           dist=repl,
                           quad=gsub('[[:digit:]]+', '', names(repl)),USA_SSHS=Cat,BASIN=basin) |>
          group_nest(USA_SSHS,quad,BASIN) |>
          left_join(mod|>select(quad,USA_SSHS,BASIN,model_asymp),,by=c("USA_SSHS","quad","BASIN"))
      }else {
        repl <- data.frame(speed=as.numeric(gsub("NE|SE|SW|NW","",names(repl))),
                           dist=repl,
                           quad=gsub('[[:digit:]]+', '', names(repl)),USA_SSHS=Cat,BASIN=basin) |>
          group_nest(quad,BASIN) |>
          left_join(mod|>filter(USA_SSHS==1)|>select(quad,BASIN,model_asymp),by=c("quad","BASIN"))
      }
      repl <- repl |>
        mutate(dist_pred=map2(model_asymp,data,predict))|>
        tidyr::unnest(c(data,dist_pred))|>
        mutate(speed=maxwind)|>
        group_nest(quad,model_asymp) |>
        mutate(dist_pred_maxwind=map2(model_asymp,data,predict))|>
        select(-(model_asymp)) |>
        tidyr::unnest(c(data,dist_pred_maxwind))|>
        mutate(dist=round(dist),
               diff=dist_pred/dist,
               dist_pred_maxwind2 = floor(dist_pred_maxwind/diff),
               name=paste0(quad,speed))|>
        select(name,dist_pred_maxwind2) |> pull(name=name)
    }else{
      ####  insert the maximum wind and its distance
      diff <- repl-maxwinddist
      repl <- replace(repl,1:length(repl),maxwinddist)
      ###  sometimes the maximum wind distance is greater than that supplied by the quadrant data
      ##  this is either due to asymmetry or sometimes seemingly to error
      ##  in these cases, make the maximum wind distance proportionally shorter than the next known highest wind distance
      repl[diff<=0] <- as.numeric(gsub('\\D+', "", names(repl[diff<=0])))/maxwind*repl[diff<=0]
      #repl[diff==0] <- repl[diff==0]-repl[diff==0]/2
      names(repl) <- gsub(gsub("NE|SE|SW|NW","",names(repl)),maxwind,names(repl))
    }
    swathpoints <- c(swathpoints,c(repl))
    #swathpoints[swathpoints<10] <-10
  }
  rot = function(a) matrix(c(round(cos(a)), round(sin(a)), -round(sin(a)), round(cos(a))), 2, 2)
  for (i in which(swathpoints>0)){
    swath <- matrix(c(0,0,
                      0,swathpoints[i]*1852,
                      swathpoints[i]*cos(80*pi/180)*1852,swathpoints[i]*sin(80*pi/180)*1852,
                      swathpoints[i]*cos(70*pi/180)*1852,swathpoints[i]*sin(70*pi/180)*1852,
                      swathpoints[i]*cos(60*pi/180)*1852,swathpoints[i]*sin(60*pi/180)*1852,
                      swathpoints[i]*cos(50*pi/180)*1852,swathpoints[i]*sin(50*pi/180)*1852,
                      swathpoints[i]*cos(40*pi/180)*1852,swathpoints[i]*sin(40*pi/180)*1852,
                      swathpoints[i]*cos(30*pi/180)*1852,swathpoints[i]*sin(30*pi/180)*1852,
                      swathpoints[i]*cos(20*pi/180)*1852,swathpoints[i]*sin(20*pi/180)*1852,
                      swathpoints[i]*cos(10*pi/180)*1852,swathpoints[i]*sin(10*pi/180)*1852,
                      swathpoints[i]*1852,0,
                      0,0),ncol=2,byrow=T)
    linestring <- swath[2:(nrow(swath)-1),]
    swath <- st_sfc(st_polygon(list(swath)),crs=custCRS)
    linestring <- st_sfc(st_linestring(linestring),crs=custCRS)
    if (grepl("SE",names(swathpoints)[i])){
      swath <- swath*rot(pi/2)
      linestring <- linestring*rot(pi/2)
    }else if (grepl("SW",names(swathpoints)[i])){
      swath <- swath*rot(pi)
      linestring <- linestring*rot(pi)
    }else if (grepl("NW",names(swathpoints)[i])){
      swath <- swath*rot(1.5*pi)
      linestring <- linestring*rot(1.5*pi)
    }
    swath <- st_as_sf(swath,crs=custCRS)
    linestring <- st_as_sf(linestring,crs=custCRS)
    swath$quad <- substr(names(swathpoints[i]),1,2)
    swath$kts <- as.numeric(gsub("NE|SE|SW|NW","",names(swathpoints)[i]))
    linestring$quad <- substr(names(swathpoints[i]),1,2)
    linestring$kts <- as.numeric(gsub("NE|SE|SW|NW","",names(swathpoints)[i]))
    swath$dist_m <- swathpoints[i]*1852
    linestring$dist_m <- swathpoints[i]*1852
    #sfCat(paste("init Iteration ", i), sep="\n")
    if (exists("swaths")){
      swaths <- rbind(swaths,swath);
      linestrings <- rbind(linestrings,linestring)
    } else{
      swaths <- swath
      linestrings <- linestring
    }
  }
  if (is.na(EYE)&Cat>0){
    eyerad <- data.frame(USA_SSHS=Cat,BASIN=basin,dist=min(swathpoints[swathpoints>0])) |>
      group_nest(USA_SSHS,BASIN) |>
      left_join(emod|>select(USA_SSHS,BASIN,model),by=c("USA_SSHS","BASIN")) |>
      mutate(pred=map2(model,data,predict))|>
      select(-model) |>
      tidyr::unnest(c(data,pred))|>
      mutate(pred=round(pred))|>pull(pred)
    if (eyerad<min(swathpoints[swathpoints>0])){
      eye = bind_rows(st_buffer(center,eyerad*1852) |>
                        st_as_sf() |>
                        mutate(quad="eye",kts=0,dist_m=eyerad*1852),
                      st_buffer(center,1) |>
                        st_as_sf() |>
                        mutate(quad="eye",kts=0,dist_m=0))
      swaths <- rbind(swaths,eye)
      linestrings <- rbind(linestrings,st_cast(eye,"LINESTRING"))
    }
  }
  if (is.na(ROCId)){
    if (Cat>0){
      ROCId <- data.frame(USA_SSHS=Cat,BASIN=basin,dist=max(swathpoints)) |>
        group_nest(USA_SSHS,BASIN) |>
        left_join(rocimod|>select(USA_SSHS,BASIN,model),by=c("USA_SSHS","BASIN")) |>
        mutate(pred=map2(model,data,predict))|>
        select(-model) |>
        tidyr::unnest(c(data,pred))|>
        mutate(pred=round(pred))|>pull(pred)
    }else{
      ROCId=max(swathpoints)+10
    }
  }
  if (ROCId<max(swathpoints)){
    ROCId <- max(swathpoints)+10
  }
  ROCI = st_buffer(center,ROCId*1852) |>
    st_as_sf()|>
    mutate(quad="ROCI",kts=0,dist_m=ROCId*1852)
  swaths <- bind_rows(swaths,ROCI)
  linestrings <- bind_rows(linestrings,st_cast(ROCI,"LINESTRING"))

  swaths <- swaths |>
    group_by(kts,quad) |>
    summarise(dist_m_mean = mean(dist_m),
              dist_m_min = min(dist_m)) |>
    ungroup() |>
    st_cast()|>
    left_join(data.frame(source=swathmodeled) |> add_rownames("quadspeed")|>rowwise()|>mutate(quad=substr(quadspeed,1,2),kts=gsub(quad,"",quadspeed)|>as.numeric()),
              by=join_by(quad,kts)) |>
    rowwise() |>
    mutate(source=if_else(!is.na(maxwinddist)&kts==maxwind,"native",if_else(is.na(source)|source==0,"modeled","native")))|>
    select(-quadspeed)
  linestrings <- suppressMessages(linestrings |>
                                    group_by(kts,quad) |>
                                    summarise(dist_m_mean = mean(dist_m),
                                              dist_m_min = min(dist_m)) |>
                                    st_cast()|>
                                    left_join(data.frame(source=swathmodeled) |> add_rownames("quadspeed")|>rowwise()|>mutate(quad=substr(quadspeed,1,2),kts=gsub(quad,"",quadspeed)|>as.numeric()),
                                              by=join_by(quad,kts)) |>
                                    rowwise() |>
                                    mutate(source=if_else(!is.na(maxwinddist)&kts==maxwind,"native",if_else(is.na(source)|source==0,"modeled","native")))|>
                                    select(-quadspeed))

  swaths$ID <- id; linestrings$ID <- id
  swaths$name <- name;linestrings$name <- name
  swaths$date <- date;linestrings$date <- date
  swaths$maxWind <- as.numeric(maxwind);linestrings$maxWind <- as.numeric(maxwind)
  swaths$minpress <- as.numeric(minpress);  linestrings$minpress <- as.numeric(minpress)
  swaths$maxpress <- as.numeric(maxpress);  linestrings$maxpress <- as.numeric(maxpress)
  swaths$centerX <- X; linestrings$centerX <- X
  swaths$centerY <- Y; linestrings$centerY <- Y
  swaths <- swaths |>
    ungroup() |>
    st_cast("POLYGON") |>
    st_transform(custCRSall) |>
    rename(geometry=x)|>
    st_set_geometry("geometry")
  linestrings <- linestrings |>
    ungroup() |>
    st_transform(custCRSall)|>
    rename(geometry=x) |>
    st_set_geometry("geometry")
  ###  include the track

  track <- tracks |>
    st_as_sf(coords=c("LON","LAT"),crs=4326,remove=FALSE)|>
    mutate(geometry_lead = lead(geometry, default = NULL)) |>
    # drop the NA row created by lagging
    slice(-n()) |>
    mutate(line = st_sfc(purrr::map2(
      .x = geometry,
      .y = geometry_lead,
      .f = ~{st_union(c(.x, .y)) |> st_cast("LINESTRING")}
    ))) |>
    st_drop_geometry()|>
    select(-geometry_lead) |>
    rename(geometry="line")|>
    st_set_geometry("geometry")|>
    st_set_crs(4326)|>
    st_transform(custCRSall)|>
    filter(ISO_TIME==track$ISO_TIME) |>
    mutate(kts=as.numeric(USA_WIND),quad="track",dist_m_mean=0,dist_m_min=0,source="native",ID=id,name=name,date=date,maxWind=as.numeric(USA_WIND),minpress=as.numeric(USA_PRES),
           maxpress=as.numeric(USA_POCI),centerX=X,centerY=Y)|>
    select(kts,quad,dist_m_mean,dist_m_min,source,ID,name,date,maxWind,minpress,maxpress,centerX,centerY)
  linestrings <- bind_rows(linestrings,track)
  cat(paste("\rBuilding Wind Extents: %",round(100*L/nrow(tracks),1)))
  return(c(swaths=list(swaths),linestrings=list(linestrings)))
}

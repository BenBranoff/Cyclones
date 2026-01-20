#' @importFrom dplyr row_number filter slice_min full_join join_by
#' @importFrom sf st_transform
#' @importFrom terra extract crs distance minmax classify values
#######  function to compare winds stored as the original IBTRaCS swaths data and the modeled high winds both stored as
########   vectors, with the interpolated raster values
compare_winds <- function(rasts,shapes=NULL,meth){
  sondedata <- get_sondes(paste0(unique(shapes$name),unique(format(shapes$date,"%Y"))),ddir="E:/OneDrive - USDA/Hurricanes/Drop Sondes")
  if (!is.data.frame(sondedata)) return(FALSE)
  sondedata$time <- as.POSIXct(sondedata$time)
  r1 <- rasts[[1]]
  if (inherits(shapes, "sf")){
  l1 <-shapes|> filter(source=="native",!quad %in% c("track","track points")) |>
    st_transform(crs(r1)) |>
    mutate(row=row_number())
  lsamps1 <- extract(r1,l1)
  lsamps <- l1 |> st_drop_geometry() |>
    full_join(lsamps1,join_by("row"=="ID")) |>
    mutate(date=unique(l1$date)) |>
    mutate(source=if_else(kts %in% c(34,50,64),"swaths","modeled")) |>
    #filter(quad!="ROCI")|>
    group_by(row,date,quad,kts,source) |>
    summarise(max=max(lyr.1,na.rm=TRUE),
              mean=mean(lyr.1,na.rm=TRUE),
              min=min(lyr.1,na.rm=TRUE)) |>
    mutate(ID=unique(shapes$ID))
  }else{;lsamps=NULL}
  rsamps <- lapply(rasts, function(x,s,m){
    sonde <- s |> mutate(difft=difftime(time(x),time,units="mins")) |> filter(abs(difft)<5) |>st_transform(crs(x))
    if (nrow(sonde)>0){
      sonde$rsamp <- extract(x,sonde)$lyr.1
      sonde$minDistWithin10kts.m <- NA
      sonde$rsource <- NA
      sonde_30ft <- which(sonde$ZW.m<10)
      if (length(sonde_30ft)>0){
        dst<- lapply(sonde_30ft,function(i,y,so) {
          #y[y<(&y>(as.numeric(so[i,]$WS.ms)*1.94-10)] <- -9999
          target_y <- classify(y, matrix(c(-Inf, as.numeric(so[i,]$WS.ms)*1.94-10, NA, as.numeric(so[i,]$WS.ms)*1.94+10, Inf, NA), ncol=3, byrow=TRUE), others=-9999)
          if (all(is.na(values(target_y)))) return(Inf)
          d <- distance(target_y)
          d <- extract(d, so[i,])[,2]
          #d <- distance(y,so[i,],rasterize=TRUE,target=-9999)
          #d <- mask(d,y)
          #d <- minmax(d)[1]
          d
        },y=x,so=sonde)
        dst <- do.call(c,dst)
        sonde$minDistWithin10kts.m[sonde_30ft] <- dst
        sonde$rsource[sonde_30ft] <- m
      }
      sonde
    }else{NULL}
  },s=sondedata,m=meth)
  rsamps <- do.call(rbind, rsamps)
  return(list(lsamps=lsamps,rsamps=rsamps))
}

#' @importFrom dplyr row_number filter slice_min full_join join_by
#' @importFrom sf st_transform
#' @importFrom terra extract crs distance minmax classify values time
#######  function to compare winds stored as the original IBTRaCS swaths data and the modeled high winds both stored as
########   vectors, with the interpolated raster values
compare_winds <- function(rasts,shapes=NULL,checksondes=TRUE){
  if (!is.null(rasts)){
    raststorm <- strsplit(names(rasts),"_")
    if (raststorm[[1]][1]=="MSW")   raststorm <- paste0(sapply(raststorm,"[[",2),sapply(raststorm,"[[",3)) |> tolower() |> unique() else
    raststorm <- paste0(sapply(raststorm,"[[",1),sapply(raststorm,"[[",2)) |> tolower() |> unique()
  }else(raststorm<-NULL)
  if (!is.null(shapes)){
    shapestorm <- paste0(unique(shapes$name),unique(format(shapes$date,"%Y")))
  }else{shapestorm<-NULL}
  if (!is.null(raststorm)&&!is.null(shapestorm)){
    storm <- ifelse(raststorm==shapestorm,raststorm,stop("Two different storms provided between raster and vector object"))
  }else if (!is.null(raststorm)){
    storm=raststorm
  }else{
    storm=shapestorm
  }
  if (checksondes) {
    for (st in raststorm){
      if (exists("sondedata", inherits=FALSE))
        sondedata=rbind(sondedata,get_sondes(st,ddir="E:/OneDrive - USDA/Hurricanes/Drop Sondes"))
      else
       sondedata <- get_sondes(st,ddir="E:/OneDrive - USDA/Hurricanes/Drop Sondes")
    }
  }
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
  rsamp_f <- function(x,s){
    m <- strsplit(names(x[[1]]),"_")[[1]][4]
    strm  = paste0(strsplit(names(x[[1]]),"_")[[1]][2],strsplit(names(x[[1]]),"_")[[1]][3]) |> tolower()
    sonde <- s |> filter(storm==strm) |> mutate(difft=difftime(time(x),time,units="mins")) |> filter(abs(difft)<5) |>st_transform(crs(x))
    if (nrow(sonde)>0){
      sonde$rsamp <- extract(x,sonde)[,2]
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
        sonde
      }
    }
  }
  if (is.list(rasts)){
    nrasts <- length(unlist(rasts))
    nlyrs <- sum(lapply(unlist(rasts),terra::nlyr)|>unlist())
    stackcount <- 0
    for (st in raststorm){
      rsts <- rasts[grep(substr(st,1,nchar(st)-4),sapply(strsplit(names(rasts),"_"),paste0,collapse="")|>tolower())] |> unlist()
      for (r in rsts){
        cat(paste0("\rSampling at drop sonde locations from ",nrasts," raster stacks and ",nlyrs," layers: %",round(100*stackcount/nrasts,1)))
        #stackvals <-  extract(r,sondedata)
        #difft <- outer(time(r), sondedata$time, FUN = "-")
        #difft <- abs(difft)<5*60
        #stackvals <- stackvals[,which(rowSums(difft)>0)+1]
        #difft <- difft[which(rowSums(difft)>0),]
        #stackvals <- stackvals[rowSums(t(difft))>0,]
        #stackvals <- stackvals[,colSums(!is.na(stackvals))>0]
        #difft <- difft[,!is.na()]

        rsamp <- lapply(r,rsamp_f,s=sondedata)
        rsamp <- lapply(Filter(Negate(is.null), rsamp),st_transform,crs=4326)
        if (exists("rsamps",inherits=FALSE)) rsamps <- rbind(rsamps,rsamp) else rsamps <- rsamp
        stackcount<-stackcount+1
      }
    }
  }else{
    rsamps <- lapply(rasts,rsampf_f,s=sondedata)
  }
  rsamps <- do.call(rbind, rsamps)
  return(list(lsamps=lsamps,rsamps=rsamps))
}

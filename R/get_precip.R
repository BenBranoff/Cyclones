#' @importFrom sf st_as_sf st_transform st_buffer st_bbox
#' @importFrom dplyr filter
#' @importFrom ecmwfr wf_request

get_precip <- function(storms,dpath=NULL){
  if (is.null(dpath)) dpath=tempdir()
  storms$ID <- ifelse(("SID" %in% names(storms)),storms$SID,storms$USA_ATCF_ID)
  storms_precip <- lapply(unique(storms$ID),function(id,y,dp){
    points <- y |> filter(ID==id) |> mutate(LAT=if_else(grepl("N",LAT),gsub("N|S| ","",LAT),paste0("-",gsub("N|S| ","",LAT))),LON=if_else(grepl("E",LON),gsub("W|E| ","",LON),paste0("-",gsub("W|E| ","",LON))))|>
      st_as_sf(coords=c('LON','LAT'),crs=4326,remove=FALSE)
    cust_CRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",mean(as.numeric(points$LON),na.rm=TRUE)," +lat_0=",mean(as.numeric(points$LAT),na..rm=TRUE))
    bounds <- points  |>
      ##  needs to be in meters to buffer the 500km radius
      st_transform(cust_CRS) |>
      ###  95% of rainfall occurs within 500km of center
      st_buffer(500000,nQuadSegs=2) |>
      mutate(IS_TIME=as.POSIXct(ISO_TIME,tz="UTC")) |>
      ###  to be able to get 48 hours before storm, subtract from mintime
      summarise(mintime=min(ISO_TIME)-48*60*60,
                maxtime=max(ISO_TIME))
    times = seq(bounds$mintime,bounds$maxtime,"hours")
    ext <- round(st_bbox(st_transform(bounds,4326)))
    ###  get the blocks of full days. These can be requested all at once
    days = unique(format(times, "%j"))
    #fulldays <- which(table(days)==24)
    #daytimes <- data.frame(day=unique(days),times=)
    precips <- lapply(days, function(d,t,e,id,dp){
      ts <- t[format(t,"%j")==d]
      request <- list(
        dataset_short_name = "reanalysis-era5-single-levels",
        product_type = "reanalysis",
        variable = "total_precipitation",#,"10m_v_component_of_wind"),
        year = unique(format(ts,"%Y")),
        month =  unique(format(ts,"%m")),
        day =  unique(format(ts,"%d")),
        time = c("00:00", "01:00", "02:00", "03:00","04:00","05:00",
                 "06:00","07:00", "08:00", "09:00", "10:00","11:00",
                 "12:00","13:00", "14:00", "15:00", "16:00","17:00",
                 "18:00","19:00","20:00", "21:00", "22:00", "23:00"),
        data_format = "grib",
        download_format = "unarchived",
        area = c(e$ymax,e$xmin,e$ymin,e$xmax),
        target = paste(id,unique(format(ts,"%Y")),unique(format(ts,"%m")),unique(format(ts,"%d")),sep="-")
      )
      if (!file.exists(paste0(dp,"/",request$target,".grib")))
      file <- wf_request(
        request  = request,  # the request
        transfer = TRUE,     # download the file
        path     = dp#paste(dp,request$target,sep="/")#paste0("C:\\Users\\BenjaminBranoff\\Downloads"       # store data in current working directory
      )
      return(rast(paste0(dp,"/",request$target,".grib")))
    },t=times,e=ext,id=id,dp=dpath)
    ###  now do the remaining days
    precips <- rast(precips)
    times <- time(precips)
    ###  find which storm segment is closest to each raster pixel
    pnts <- as.data.frame(precips[[1]],xy=TRUE)
    pnts <- pnts%>%st_as_sf(coords=c("x","y"),crs=crs(precips[[1]]))%>%st_transform(cust_CRS )
    ### rains start about 500km out, so find the ctrack that first gets within 500km of a point
    ###  then find the last track that is within 500km of a point, those are the start and end points for the storm precip.
    rains <- points %>%
      st_transform(cust_CRS) %>% st_buffer(500000)
    rains <- st_covered_by(pnts,rains)
    rainstart <- lapply(rains,function(x) ifelse(length(x)==0,NA,min(x)))
    rainsend <- lapply(rains,function(x) ifelse(length(x)==0,NA,max(x)))
    pnts$start <- unlist(rainstart)
    pnts$end <- unlist(rainsend)
    pnts <- st_transform(pnts,crs(precips))
    start <- rasterize(pnts,precips[[1]],field="start")
    end <- rasterize(pnts,precips[[1]],field="end")
    #pnts$closest <- alltracks%>%filter(SID==storm)%>%
    #  st_as_sf(coords=c("LON","LAT"),crs=4326) %>% st_nearest_feature(pnts,.)
    #closest <- rasterize(pnts,precips[[1]],field="closest")
    ###  which of the storm tracks times matches the precipitation times
    matchtimes <-  match(points%>% pull(ISO_TIME),times)
    while(sum(is.na(matchtimes))>0){
      matchtimes[is.na(matchtimes)] <- lag(matchtimes)[is.na(matchtimes)]
    }
    ###  use the closest storm segment to get the previous week of precipitation files
    ##  some 3hr windows are missing, so get the actual time difference between each file
    # difftimes <- difftime(lead(times),times,"hour")
    #difftimes[is.na(difftimes)] <- 3
    start_pre7 <- app(start,function(x,y) y[x]-24*7/3,y=matchtimes)
    start_pre2 <- app(start,function(x,y) y[x]-24*2/3,y=matchtimes)
    start_storm <- end_pre <- app(start,function(x,y) y[x],y=matchtimes)
    end_storm <- app(end,function(x,y) y[x],y=matchtimes)
    ### to get hourly rate, multiply by the number of hours between each
    #precips2 <- precips/as.numeric(difftimes)
    precip_presum7 <- rapp(precips, start_pre7, end_pre, sum,na.rm=TRUE)
    precip_presum2 <- rapp(precips, start_pre2, end_pre, sum,na.rm=TRUE)
    precip_stormsum <- rapp(precips, start_storm,end_storm,sum,na.rm=TRUE)
    precip_max= rapp(precips,start_storm,end_storm,max,na.rm=TRUE)
    precip_mean= rapp(precips,start_storm, end_storm, mean,na.rm=TRUE)
    precips <- c(precip_presum7,precip_presum2,precip_stormsum,precip_max,precip_mean)
    names(precips) = paste0(id,c("sum_prestorm7d","sum_prestorm2d","sum_storm","max","mean"))
    #precips <- mask(precips,alltracks%>%filter(SID==storm)%>%st_as_sf(coords=c("LON","LAT"),crs=4326)%>%st_transform(crs(precips))%>%st_buffer(500000)%>%st_union%>%st_as_sf())
    precips
  } ,y=storms)
  storms_precip
}


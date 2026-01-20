#' @importFrom sf st_length st_line_sample st_crs st_geometry
#' @keywords internal
interp_track <- function(track,t_res){
  ###  this function interpolates points along a line segment according to a time interval, not a spacing interval
  ###  this assumes the speed of the storm is constant during the segment duration (typically 3 hours)
  ###  get the end points of the track segments
  points <- track |> filter(quad=="track points")
  track <- track |> filter(quad=="track")
  #points <- track|>
  #  st_cast("POINT")|>
  #  ### take only the first end point of each line
  #  group_by(date)|>
  #  slice(1)|>
  #  select(date) |>
  #  rename(ISO_TIME=date)|>
  #  ungroup()
  #if (length(track$STORM_SPEED)==0){
  #  track <- track |>
  #    ###  meters per second to knots
  #    mutate(speed=1.94384449*as.numeric(st_length(track))/as.numeric(difftime(date,lag(date),unit="secs")),
  #           STORM_SPEED=if_else(is.na(speed),lead(speed),speed))
  #}
  ####  for each segment, add new points along the length of the segment
  ###  points should be spaced equally apart according to location, not time
  for (n in 1:(nrow(track))){
    ## create custom CRS centered on segment to minimize geometry calculation errors
    custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",track$centerX[n]," +lat_0=",track$centerY[n])
    t1 <- track$date[n]
    if (n==nrow(track)){
    t2 <- points$date[n+1]
    }else{
    t2 <- track$date[n+1]}
    ### if this is the last segment of the track,we dont know its duration because there is no next track with a timestamp
    ### need to calculate its duration, which we can do with the given speed and length

      #if (track$STORM_SPEED[nrow(track)]==0){
      #  speed=0.1
      #}else{
      #  speed=as.numeric(track$STORM_SPEED[nrow(track)])
      #}
      ####  the length is in meters but the speed is in knots
      ###  one knot is 0.514444 meters per second
      ###  this calculated dt is in hours
    #   dt <- round(as.numeric(st_length(track[nrow(track),]))/(speed*0.5144444),1)/(60)
    #   t2 <- t1+(1+round(dt))*60*60
    #   ####  create the final end point of the storm
    #   end_point <- track[n,] |>
    #     st_transform(custCRS) |>
    #     st_line_sample(type="regular",sample=ints) |>
    #     st_as_sf(crs=custCRS) |>
    #     st_transform(st_crs(track))|>
    #     st_cast("POINT")|>
    #     mutate(ISO_TIME =t2)
    #   sf::st_geometry(end_point)="geometry"
    # }else{
      ###  get the time difference between the two adjacent segments
      dt = as.numeric(difftime(t2,t1,units="mins"))
   # }
    ###  only necessary to interpolate if they are more than the desired temporal resolution
    if (dt>t_res){
      ###  establish the times we need in between the two known locations
      ###  in this case, in 3 minute intervals
      int_times <- c(t1,
                     seq.POSIXt(as.POSIXct(ceiling(as.numeric(t1)/(t_res*60))*t_res*60,origin = "1970-01-01",tz="GMT"),#lubridate::ceiling_date(t1,"3 mins"),
                                as.POSIXct(floor(as.numeric(t2)),tz="GMT"),by=t_res*60),#lubridate::floor_date(t2),by="3 mins"),
                     t2)
      ###  remove duplicates often created at extremes of the interval
      int_times <- int_times[!duplicated(int_times)]
      ###  get the time difference between the interpolated times
      dt <- as.numeric(difftime(int_times,dplyr::lag(int_times),units="hours"))
      ##  if there are only two points, just take the end point
      if (length(dt)<3){
        dt <- dt[2]
        ints <- dt
        ###  otherwise, take everything but the end point
        ##  the intervals are the cumulative proportion of the segment as calculated by the
        ##  dt time inervals relative to the duration of the entire segment
      }else{
        dt <- dt[2:(length(dt)-1)]
        ints <- cumsum(dt[!is.na(dt)])/as.numeric(difftime(t2,t1,units="hours"))
      }
      ###  now sample along the track segment according to the intervals calculated
      ###  transform back to original crs
      t  <- track[n,] |>
        st_transform(custCRS) |>
        st_line_sample(type="regular",sample=ints) |>
        st_as_sf(crs=custCRS) |>
        st_transform(st_crs(track))|>
        st_cast("POINT")|>
        mutate(date = int_times[2:(length(int_times)-1)])#track$ISO_TIME[n]+60*15*seq(1,length(dt)))
      sf::st_geometry(t)="geometry"
      if (n==1){newpoints <- t}else{newpoints <- rbind(newpoints,t)}
    }
  }
  ##  add all together: segment end points, interpolated points and final end point
  points <- rbind(points |> select(date),newpoints)|>
    dplyr::arrange(date)
  ###  find the outer bounds of the storm
  points
}

#' @importFrom dplyr select mutate rename filter bind_rows if_else across group_by ungroup summarise nest_by reframe distinct group_nest left_join
#' @importFrom tidyr pivot_longer unnest
#' @importFrom purrr map map2
#' @importFrom stats dist lm nls predict setNames
build_models <- function(tracks=NULL){
  if (is.null(tracks))
  tracks <- get_storms(source="ncei",ib_filt=NULL)
  ##  collect wind extent observations from the data
  ##  These are supplied as the extent of 34, 50 and 64 knot sustained wind speeds within each of the 4 cardinal direction quadrants of the storm (NE,SE,SW,NW)
  ##  also, for some storms the maximum wind AND the extent of maximum winds is provided
  quaddist <- bind_rows(
    tracks |>
      ####  these columns describe the extent of 34,50 and 64 knot winds in each of the storms quadrant
      ####  at different timesteps
      select(BASIN,USA_R34_NE,USA_R34_SE,USA_R34_NW,USA_R34_SW,
             USA_R50_NE,USA_R50_SE,USA_R50_NW,USA_R50_SW,
             USA_R64_NE,USA_R64_SE,USA_R64_NW,USA_R64_SW,
             ###  these are identifiers and more variables important for additional models
             SID,ISO_TIME,USA_SSHS,USA_EYE,USA_RMW,USA_WIND,USA_PRES,USA_ROCI,USA_POCI,USA_PRES,BASIN) |>
      pivot_longer(cols=c(USA_R34_NE,USA_R34_SE,USA_R34_NW,USA_R34_SW,
                                 USA_R50_NE,USA_R50_SE,USA_R50_NW,USA_R50_SW,
                                 USA_R64_NE,USA_R64_SE,USA_R64_NW,USA_R64_SW),
                          names_to =c("speed","quad"),
                          names_pattern = "USA_R(.*)_(.*)",values_to="dist") |>
      mutate(speed=as.numeric(speed),
             source="swaths"),
    tracks |>
      select(BASIN,SID,ISO_TIME,USA_SSHS,USA_EYE,USA_RMW,USA_WIND,USA_PRES,USA_ROCI,USA_POCI,USA_PRES,BASIN) |>
      ###  there are also two columns that describe the maximum wind speed and its maximum distance in the storm
      ###  this is not provided for each quadrant and so is assumed to be circular and symmetric
      mutate(quad=list(c("NE","SE","SW","NW")),
             speed=as.numeric(USA_WIND),
             source="RMW") |>
      select(-USA_WIND) |>
      unnest(quad) |>
      rename(dist=USA_RMW)) |>
    mutate(across(c(dist, speed,USA_EYE,USA_ROCI,USA_POCI,USA_PRES,USA_SSHS), ~as.numeric(.)),
           USA_SSHS_lab = paste("Category ",USA_SSHS),
           BASIN=if_else(is.na(BASIN),"NA",BASIN)) |>
    filter(dist>0)
  ############################################################
  #######   create models needed later
  ######    We need these to reconstruct maximum wind extents from given data
  ########################################################
  ###  a non-linear model of the wind speed distance as a function of wind speed and storm size by quadrant
  ###  in general, the SW quadrant has higher winds closer to the center, and the NE quadrant has the same wind speeds farther out
  ###  storm size also seems to be important
  ###  Asymptotic models by storm size and quadrant
  modquad <- lm(data=quaddist|>filter(USA_SSHS>=1),dist~log(speed):(quad*BASIN))
  modquad_asymp <- nls(dist ~ SSasymp(speed, Asym,R0, lrc), data = quaddist|> filter(USA_SSHS>=1),control = list(maxiter = 500))
  ###  try the model grouped by storm size (saffir simpson)
  ###  create a new dataset to predict from the model, with all possibilities of wind speed
  ##  create new dataset with all wind speeds present for each storm size and quadrant
  quaddist_new <- quaddist |>
    filter(USA_SSHS>=1,BASIN!="SA")|>
    distinct(USA_SSHS,quad,BASIN)|>
    group_by(USA_SSHS,quad,BASIN) |>
    reframe(speed=1:185)|>
    nest_by(USA_SSHS,quad,BASIN)
  ##  build models on the original data, nested by both size, quadrant, and BASIN
  modquads <- quaddist|> filter(USA_SSHS>=1) |>
    group_nest(USA_SSHS,quad,BASIN) |>
    ## south atlantic has too few storms
    filter(BASIN!="SA")|>
    mutate(data2=quaddist_new$data,  ###  predict fitted models onto the new dataset with the full range of windspeeds
           model_asymp = map(data,~nls(dist ~ SSasymp(speed, Asym,R0, lrc), data = .,control = list(maxiter = 500))),
           model_lm = map(data,~lm(dist ~log(speed), data = .)),
           pred_asymp=map2(model_asymp,data2,predict),
           pred_lm=map2(model_lm,data2,predict))
  ### take the predicted and original values to plot with
  modquad_df <- modquads |>
    select(-c(data,model_asymp,model_lm)) |>
    unnest(c(data2,pred_asymp,pred_lm)) |>
    left_join(modquads|> unnest(data)|>select(USA_SSHS,quad,BASIN,speed,dist),by=c("USA_SSHS","quad","BASIN","speed")) |>
    distinct(USA_SSHS,quad,BASIN,speed,dist,.keep_all=T)
  ####  the non-linear model is a better fit for all combinations of quad and saffir simpson
  #####
  ##  also need to model the eyewall radius
  ####
  eyemod <- quaddist |>
    group_by(SID,ISO_TIME) |>
    filter(dist==min(dist,na.rm=T),USA_EYE!=" ")|>
    mutate(USA_EYE=as.numeric(USA_EYE)/2) |>
    ungroup() |>
    group_nest(USA_SSHS,BASIN) |>
    mutate(model=map(data,~lm(USA_EYE~0 + dist,data=.)),
           pred=map(model, predict))
  ##  and storms size
  ####
  ###  for larger storms, the size can be reasonably predicted from the largest known wind swath
  ROCImod <- quaddist |>
    group_by(SID,ISO_TIME) |>
    filter(dist==max(dist,na.rm=T),!is.na(USA_ROCI)) |>
    ungroup() |>
    group_nest(USA_SSHS,BASIN) |>
    mutate(model=map(data,~lm(USA_ROCI~dist,data=.)),
           pred=map(model, predict))
  ###  also need to predict POCI when it is absent
  POCImod <- quaddist |>
    filter(!is.na(USA_POCI),!is.na(USA_PRES),(USA_POCI-USA_PRES)>0) |>
    group_by(SID,ISO_TIME) |>
    mutate(USA_WIND=max(speed),
           Pressdif = USA_POCI-USA_PRES) |>
    ungroup()|>
    group_nest(USA_SSHS,BASIN) |>
    mutate(model=map(data,~lm(Pressdif~USA_WIND,data=.)),
           pred=map(model, predict))
  ## also need the average minimum pressure by month, basin, and storm size for rare cases when missing
  minpress <- quaddist |> mutate(MONTH=format.Date(ISO_TIME,"%m")) |>group_by(BASIN,USA_SSHS,MONTH) |>summarise(USA_PRES=mean(USA_PRES,na.rm=TRUE))
  mods <- list(modquads=modquads|> select(-c(data,pred_asymp,pred_lm)),
               emod=eyemod|> select(-c(data,pred)),
               rocimod=ROCImod|> select(-c(data,pred)),
               pocimod=POCImod|> select(-c(data,pred)),
               minpress = minpress)

 mods
}

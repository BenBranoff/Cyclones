#' @importFrom terra extract
#######  function to compare winds stored as the original IBTRaCS swaths data and the modeled high winds both stored as
########   vectors, with the interpolated raster values
compare_winds <- function(rasts,shapes){
  l1 <-shapes[[1]] %>%
    st_transform(crs(rasts[[1]]$vel)) %>%
    mutate(row=row_number())
  l2 <- shapes[[2]] %>%
    st_transform(crs(rasts[[2]]$vel)) %>%
    mutate(row=row_number())
  lsamps1 <- extract(rasts[[1]]$vel,l1)
  lsamps2 <- extract(rasts[[2]]$vel,l2)
  lsamps <- bind_rows(
    l1 %>% st_drop_geometry() %>%
      full_join(lsamps1,join_by("row"=="ID")) %>%
      mutate(date=unique(l1$date)),
    l2 %>% st_drop_geometry() %>%
      full_join(lsamps2,join_by("row"=="ID"))%>%
      mutate(date=unique(l2$date)))%>%
    mutate(source=if_else(kts %in% c(34,50,64),"swaths","modeled")) %>%
    #filter(quad!="ROCI")%>%
    group_by(row,date,quad,kts,source) %>%
    summarise(max=max(lyr.1,na.rm=TRUE),
              mean=mean(lyr.1,na.rm=TRUE),
              min=min(lyr.1,na.rm=TRUE)) %>%
    mutate(ID=unique(shapes[[1]]$ID))
  return(lsamps)
}

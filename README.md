Cyclones Introduction
================
Ben Branoff
2025-03-18

R utilities for modeling and analyzing meteorological cyclones. Fit
models to observations of storm characteristics and use these models to
predict and interpolate wind and pressure fields across a storm’s area.
\## Build Wind Models
<details>
<summary>
\## Read in and parse previous cyclone data
</summary>

Start with the IBTrACS data (updated data available at:
<https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/>)

``` r
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

###  Load in the North Atlantic data only:
IBTrACS <- read.csv(paste0(dirname(getwd()),"./data/ibtracs.NA.list.v04r01.csv"))
####  limit to storms in or after 2000 for brevity
IBTrACS <- IBTrACS %>%
  filter(SEASON>=2000)
```

Determine the extent of maximum winds and other wind extent observations
from the data

``` r
quaddist <- bind_rows(
  IBTrACS %>%  
    ####  these columns describe the extent of 34,50 and 64 knot winds in each of the storms quadrant
    ####  at different timesteps 
    select(USA_R34_NE,USA_R34_SE,USA_R34_NW,USA_R34_SW,  
           USA_R50_NE,USA_R50_SE,USA_R50_NW,USA_R50_SW,  
           USA_R64_NE,USA_R64_SE,USA_R34_NW,USA_R64_SW,
           ###  these are identifiers and more variables important for additional models
           SID,ISO_TIME,USA_SSHS,USA_EYE,USA_RMW,USA_WIND,USA_PRES,USA_ROCI,USA_POCI,USA_PRES,BASIN) %>%
    tidyr::pivot_longer(cols=c(USA_R34_NE,USA_R34_SE,USA_R34_NW,USA_R34_SW,  
                               USA_R50_NE,USA_R50_SE,USA_R50_NW,USA_R50_SW,  
                               USA_R64_NE,USA_R64_SE,USA_R34_NW,USA_R64_SW),  
                        names_to =c("speed","quad"),  
                        names_pattern = "USA_R(.*)_(.*)",values_to="dist") %>%  
    mutate(speed=as.numeric(speed),  
           source="swaths"),  
  IBTrACS %>%  
    select(SID,ISO_TIME,USA_SSHS,USA_EYE,USA_RMW,USA_WIND,USA_PRES,USA_ROCI,USA_POCI,USA_PRES,BASIN) %>%
    ###  there are also two columns that describe the maximum wind speed and its maximum distance in the storm
    ###  this is not provided for each quadrant and so is assumed to be circular and symmetric
    mutate(quad=list(c("NE","SE","SW","NW")),  
           speed=as.numeric(USA_WIND),  
           source="RMW")%>%  
    dplyr::select(-USA_WIND) %>%  
    tidyr::unnest(quad) %>%  
    rename(dist=USA_RMW)) %>%  
    mutate(across(c(dist, speed,USA_EYE,USA_ROCI,USA_POCI,USA_PRES), ~as.numeric(.)),
           USA_SSHS_lab=paste0("Cat. ",USA_SSHS))
```

</details>

## Fit models from previous hurricane data

We need these to reconstruct wind extents for any given storm with
limited information. Here are non-linear asymptotic models as well as
log-tansformed linear models of the wind speed distance as a function of
wind speed and storm size by quadrant:

``` r
###  Asymptotic models by storm size and quadrant
ggplot(quaddist%>%filter(!is.na(dist)),aes(x=speed,y=dist,col=factor(quad)))+
  geom_point(size=0.5)+
  geom_point(data=quaddist%>%group_by(quad,speed,USA_SSHS) %>%summarise(dist=mean(dist,na.rm=T)))+
  stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc),se=FALSE)+
  facet_wrap(~USA_SSHS)+
  labs(col="Quadrant")+
  ggtitle("North Atlantic Tropical Cyclone Wind Speed Distance by Storm Saffir-Simpson Scale\nAsymptotic models")+
  xlab("Wind Speed (kts)")+ylab("Max. distance from center (nmi)")
```

![](README_files/figure-gfm/plot%20models-1.png)<!-- -->

``` r
####  close up of Category 5 storms
ggplot(quaddist%>%filter(!is.na(dist),USA_SSHS==5),aes(x=speed,y=dist,col=factor(quad)))+
  geom_point(size=0.5)+
  geom_point(data=quaddist%>%group_by(quad,speed,USA_SSHS) %>%summarise(dist=mean(dist,na.rm=T)),size=1)+
  stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc),se=FALSE)+
  facet_wrap(~USA_SSHS_lab)+
  labs(col="Quadrant")+
  ggtitle("North Atlantic Tropical Cyclone Wind Speed Distance by Storm Saffir-Simpson Scale\nAsymptotic models")+
  xlab("Wind Speed (kts)")+ylab("Max. distance from center (nmi)")
```

![](README_files/figure-gfm/plot%20models-2.png)<!-- -->

``` r
###  a linear model of the wind speed distance as a function of wind speed and storm size
###  fit doesnt seem as good as the non-linear model
ggplot(quaddist,aes(x=log(speed),y=dist,col=quad))+
geom_point()+
  geom_point(data=quaddist%>%
               filter(!is.na(dist))%>%
               group_by(quad,USA_SSHS,speed) %>%
               summarise(dist=mean(dist,na.rm=T))) +
  stat_smooth(method="lm", formula=y~log(x),se=TRUE)+
  facet_wrap(~USA_SSHS)+
  labs(col="Quadrant")+
  ggtitle("North Atlantic Tropical Cyclone Wind Speed Distance by Storm Saffir-Simpson Scale\nLinear models")+
  xlab("Wind Speed (kts)")+ylab("Max. distance from center (nmi)")
```

![](README_files/figure-gfm/plot%20models-3.png)<!-- --> In general, the
SW quadrant has higher winds closer to the center, and the NE quadrant
has the same wind speeds farther out. Storm size also seems to be
important. The linear model does not seem to fit as well as the
asymptotic models. In either case, models generally only compute for
true hurricanes/typhoons with Saffir-Simpson scale larger than 0. We
cant easily include an interaction term for an Asymp model, but that is
approached in further steps.

## Compare models

``` r
modquad <- lm(data=quaddist%>%filter(USA_SSHS>=1),dist~log(speed):quad)  
modquad_asymp <- nls(dist ~ SSasymp(speed, Asym,R0, lrc), data = quaddist%>% filter(USA_SSHS>=1),control = list(maxiter = 500))
summary(modquad_asymp)
```

    ## 
    ## Formula: dist ~ SSasymp(speed, Asym, R0, lrc)
    ## 
    ## Parameters:
    ##       Estimate Std. Error t value Pr(>|t|)    
    ## Asym  13.77095    0.42665   32.28   <2e-16 ***
    ## R0   777.11506   14.80396   52.49   <2e-16 ***
    ## lrc   -2.85674    0.01092 -261.68   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 39.09 on 73335 degrees of freedom
    ## 
    ## Number of iterations to convergence: 5 
    ## Achieved convergence tolerance: 6.984e-07
    ##   (10212 observations deleted due to missingness)

``` r
summary(modquad)
```

    ## 
    ## Call:
    ## lm(formula = dist ~ log(speed):quad, data = quaddist %>% filter(USA_SSHS >= 
    ##     1))
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -98.02 -27.84  -7.47  19.95 364.66 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        463.4697     1.6356   283.4   <2e-16 ***
    ## log(speed):quadNE  -98.7214     0.4117  -239.8   <2e-16 ***
    ## log(speed):quadNW -100.2221     0.4181  -239.7   <2e-16 ***
    ## log(speed):quadSE  -99.9373     0.4117  -242.7   <2e-16 ***
    ## log(speed):quadSW -103.6325     0.4118  -251.7   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 41.11 on 73333 degrees of freedom
    ##   (10212 observations deleted due to missingness)
    ## Multiple R-squared:  0.4646, Adjusted R-squared:  0.4646 
    ## F-statistic: 1.591e+04 on 4 and 73333 DF,  p-value: < 2.2e-16

``` r
AIC(modquad,modquad_asymp)
```

    ##               df      AIC
    ## modquad        6 753213.5
    ## modquad_asymp  4 745809.1

``` r
BIC(modquad,modquad_asymp)
```

    ##               df      BIC
    ## modquad        6 753268.7
    ## modquad_asymp  4 745846.0

By AIC and BIC, the asymptotic model is marginally better, not enough
for a true determination of fitness. Try new models grouped by storm
size (saffir simpson and or quad, rather than a single model with these
as covariables (in the linear model’s case).

``` r
##  create new dataset with all wind speeds present for each storm size and quadrant
quaddist_new <- quaddist %>%
  filter(USA_SSHS>=1)%>%
  distinct(USA_SSHS,quad)%>%
  group_by(USA_SSHS,quad) %>%
  reframe(speed=1:175)%>%
  nest_by(USA_SSHS,quad) 

##  build models on the original data, nested by both size and quadrant
modquads <- quaddist%>% filter(USA_SSHS>=1) %>%
  group_nest(USA_SSHS,quad) %>%
  mutate(data2=quaddist_new$data,  ###  predict fitted models onto the new dataset with the full range of windspeeds
         model_asymp = map(data,~nls(dist ~ SSasymp(speed, Asym,R0, lrc), data = .,control = list(maxiter = 500))),
         model_lm = map(data,~lm(dist ~log(speed), data = .)),
         pred_asymp=map2(model_asymp,data2,predict),
         pred_lm=map2(model_lm,data2,predict)) 
### take the predicted and original values to plot with
modquad_df <- modquads %>%
  select(-c(data,model_asymp,model_lm)) %>%
  tidyr::unnest(c(data2,pred_asymp,pred_lm)) %>%
  left_join(modquads%>% tidyr::unnest(data)%>%select(USA_SSHS,quad,speed,dist),by=c("USA_SSHS","quad","speed")) %>%
  distinct(USA_SSHS,quad,speed,dist,.keep_all=T)

ggplot(modquad_df %>% group_by(USA_SSHS,quad,speed) %>%
         summarise(dist=mean(dist,na.rm=T),
                   speed=mean(speed),
                   fit_lm=mean(pred_lm),
                   fit_asymp=mean(pred_asymp)),
       aes(x=speed,y=dist,col=factor(quad)))+
  geom_point()+
  facet_wrap(~USA_SSHS)+
  geom_line(aes(x=speed,y=fit_asymp,lty="asympt."))+
  geom_line(aes(x=speed,y=fit_lm,lty="linear"))+
  scale_linetype_manual(values=c(1,2))+
  labs(lty="model type",col="Quadrant")+
  ggtitle("North Atlantic Tropical Cyclone Wind Speed Distance by Storm Saffir-Simpson Scale\nGrouped models")+
  xlab("Wind Speed (kts)")+ylab("Max. distance from center (nmi)")
```

![](README_files/figure-gfm/grouped%20models-1.png)<!-- -->

``` r
##compare the models
bind_rows(map2(modquads$model_asymp,modquads$model_lm,~AIC(.x,.y)) %>%
                        bind_rows() %>%
                        rename(score=AIC) %>%
                        mutate(USA_SSHS=rep(modquads$USA_SSHS,each=2),
                               quad=rep(modquads$quad,each=2),
                               model=rep(c("asymp","lm"),20),
                               type="AIC"),
                      map2(modquads$model_asymp,modquads$model_lm,~BIC(.x,.y)) %>%
                        bind_rows() %>%
                        rename(score=BIC) %>%
                        mutate(USA_SSHS=rep(modquads$USA_SSHS,each=2),
                               quad=rep(modquads$quad,each=2),
                               model=rep(c("asymp","lm"),20),
                               type="BIC"))%>%
  group_by(USA_SSHS,quad) %>%
  mutate(bestAIC=model[score==min(score[type=="AIC"])],
         bestBIC=model[score==min(score[type=="BIC"])])
```

    ## # A tibble: 80 × 8
    ## # Groups:   USA_SSHS, quad [20]
    ##       df  score USA_SSHS quad  model type  bestAIC bestBIC
    ##    <dbl>  <dbl>    <int> <chr> <chr> <chr> <chr>   <chr>  
    ##  1     4 98266.        1 NE    asymp AIC   asymp   asymp  
    ##  2     3 98646.        1 NE    lm    AIC   asymp   asymp  
    ##  3     4 74989.        1 NW    asymp AIC   asymp   asymp  
    ##  4     3 75260.        1 NW    lm    AIC   asymp   asymp  
    ##  5     4 97412.        1 SE    asymp AIC   asymp   asymp  
    ##  6     3 97743.        1 SE    lm    AIC   asymp   asymp  
    ##  7     4 96399.        1 SW    asymp AIC   asymp   asymp  
    ##  8     3 96734.        1 SW    lm    AIC   asymp   asymp  
    ##  9     4 40736.        2 NE    asymp AIC   asymp   asymp  
    ## 10     3 41063.        2 NE    lm    AIC   asymp   asymp  
    ## # ℹ 70 more rows

The asymptotic model is a better fit for most combinations of quadrant
and saffir simpson scale. More importantly, the asymptotic model never
drops below zero, which is important for predicting high winds near the
eye of the storm.

Also need to model the eyewall radius. For these, its important to set
the intercept to 0. We know that the eyewall radius can’t be negative
and it cant be more than the extent of maximum winds, and setting the
intercept to 0 assures that.

``` r
eyemod <- quaddist %>% 
  group_by(SID,ISO_TIME) %>%
  filter(dist==min(dist,na.rm=T),USA_EYE!=" ")%>%
  mutate(USA_EYE=as.numeric(USA_EYE)/2) %>%
  ungroup() %>%
  group_nest(USA_SSHS) %>%
  mutate(model=map(data,~lm(USA_EYE~0 + dist,data=.)),
         pred=map(model, predict)) 


ggplot(data=eyemod %>% tidyr::unnest(c(data,pred)),
       aes(x=dist,y=USA_EYE))+
  geom_point()+
  facet_wrap(~USA_SSHS)+
  stat_smooth(method="lm",formula="y~x+0")+
  geom_line(aes(y=pred),lty=2)+
  geom_abline(slope=1,intercept=0)+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  ggtitle("North Atlantic Tropical Cyclone Eye Wall Distance and Distance of Maximum Wind by Storm Saffir-Simpson Scale\nGrouped models")+
  xlab("Minimum Distance of Maximum Wind (nmi)")+ylab("Eyewall Radius (nmi)")
```

![](README_files/figure-gfm/eyemodels-1.png)<!-- -->

These models look decent for “large” storms (2-5) but begin to weaken
considerably for “small” storms (0 & 1). We should use them cautiously
and scrutinize any predictions made.

Finally, we need models for the overall storm size or its outer radius.
This is important for constraining the interpolations made later. The
storm size is coded in IBTrACS as the radius of last closed isobar, or
ROCI. Its pressure is codes as the pressure at the last closed isobar,
or POCI.

``` r
ROCImod <- quaddist %>% 
  group_by(SID,ISO_TIME) %>%
  filter(dist==max(dist,na.rm=T),!is.na(USA_ROCI)) %>%
  ungroup() %>%
  group_nest(USA_SSHS) %>%
  mutate(model=map(data,~lm(USA_ROCI~dist,data=.)),
         pred=map(model, predict))
ggplot(data=ROCImod %>% tidyr::unnest(c(data,pred)),
       aes(x=dist,y=USA_ROCI))+
  geom_point()+
  facet_wrap(~USA_SSHS)+
  stat_smooth(method="lm",formula="y~x")+
  geom_line(aes(y=pred),lty=2)+
  geom_abline(slope=1,intercept=0)+
  ggtitle("North Atlantic Tropical Cyclone Outer Isobar Distance and Maximum Swath Distance\nGrouped models")+
  xlab("Maximum Distance of Minimum Wind (nmi)")+ylab("Radius of Last Closed Isobar (nmi)")
```

![](README_files/figure-gfm/rocimodels-1.png)<!-- -->

For larger storms, the outer radius/diameter can be reasonably predicted
from the largest known wind swath. We also need the minimum pressure of
the system. This is important for predicting the pressure at any point
in the storm, which is imortant for approximating the density of air,
which gets us to the force or that air blowing at a particulat speed,
which finally gives us the energy dissipated by that mass of air on a
surface.

``` r
POCImod <- quaddist %>% 
  filter(!is.na(USA_POCI),!is.na(USA_PRES),(USA_POCI-USA_PRES)>0) %>%
  group_by(SID,ISO_TIME) %>% 
  mutate(USA_WIND=max(speed),
         Pressdif = USA_POCI-USA_PRES) %>%
  ungroup()%>%
  group_nest(USA_SSHS) %>%
  mutate(model=map(data,~lm(Pressdif~USA_WIND,data=.)),
         pred=map(model, predict))
```

``` r
ggplot(POCImod %>% tidyr::unnest(c(data,pred)),
       aes(x=USA_WIND,y=Pressdif))+
  geom_point()+
  facet_wrap(~USA_SSHS)+
  stat_smooth(method="lm",formula=y~x)+
  geom_line(aes(y=pred),lty=2)+
  ggtitle("North Atlantic Tropical Cyclone Pressure Difference and Maximum Wind Speed\nGrouped models")+
  xlab("Maximum Wind Speed (kts)")+ylab("Storm Pressure Difference (mbar)")
```

![](README_files/figure-gfm/pocimodplot-1.png)<!-- -->

With these models, we can reconstruct a cyclone’s wind and pressure
fields from limited information. This allows us to predict velocity,
pressure, and power dissipation at any point in a storm based on
observations from previous storms. In most cases, these models are only
reasonable for true hurricanes/typhoons with a Saffir-Simpson rating
greater than 0. Save the models to use in later steps.

    saveRDS(modquads%>% select(-c(data,pred)),"modquads.rds")
    saveRDS(eyemod%>% select(-c(data,pred)),"eyemod.rds")
    saveRDS(ROCImod%>% select(-c(data,pred)),"ROCI.rds")
    saveRDS(POCImod%>% select(-c(data,pred)),"POCImod.rds")
    minpress <- quaddist %>% mutate(MONTH=format.Date(ISO_TIME,"%m")) %>%group_by(BASIN,USA_SSHS,MONTH) %>%summarise(USA_PRES=mean(USA_PRES,na.rm=TRUE))
    write.csv(minpress,"minpress.csv",row.names=FALSE)
    ###  note, when reading this data set back in, make sure the NA basin for North America is not set to missing (NA)
    ##  minpress <- minpress %>% mutate(BASIN=if_else(is.na(BASIN),"NA",BASIN))
    
## Reconstruct wind and pressure fields for individual storms

With the models from the above steps, we can now predict the extent of different wind speeds, including the maximum winds, as well as the radius of the eye wall, and the maximum extent of a storm, for any previous storm that does not include this information. 

Hurricane Camille made landfall in the Gulf Coast state of Mississippi in August of 1969. It was a devastating storm in many ways, but the IBTRACs data does not include information on the extent of any of the wind speeds. It includes only the maximum wind speed and one observation with the extent. To properly assess the intensity of the storm at any point along its path, we need to first predict where the different wind speeds (34, 50, and 64 knots) occurred. Here, we predict where the different winds occurred based on the minimum pressure and the maximum wind speed, and then project the location of the various wind speeds (34, 50, and 64 knots) onto arcs described by the predicted radius and quadrant. We do this iteratively for each time step. 
``` r
###  we will need some additional libraries
###  if we want to use parallel processing, we will need to load the snowfall package as well
library(sf)
library(dplyr)
library(purrr)
library(tibble)
#library(snowfall)

###  grab the tabular data for Camille
IBTrACS <- read.csv("ibtracs.ALL.list.v04r01.csv") %>%
  filter(SEASON==1969,NAME=="CAMILLE")%>%
  ##  make sure the timestamp is appropriately formatted
  mutate(ISO_TIME = as.POSIXct(ISO_TIME,tz="GMT"),
  MONTH=format.Date(ISO_TIME,"%m"),
  BASIN=if_else(is.na(BASIN),"NA",BASIN))
  
##source the functions
source("Swath Maker helpers.R")

##  use the 'make_swaths' function interatively over each timestep to reconstruct the spatial extent of all wind speeds
linesswaths <- lapply(1:nrow(IBTrACS),make_swaths,mod=modquads,emod=eyemod,rocimod=ROCImod,pocimod=POCImod,tracks=IBTrACS,minpresss=minpress)

##  the parallel version. Only use if confident in machine's ability to do so
##  choose an appropriate number of cpus (roughly %50-75 of available cpus)
#sfInit(parallel=TRUE, cpus=15)
#load the appropiate libraries
#sfLibrary(sf)
#sfLibrary(dplyr)
#sfLibrary(purrr)
#sfLibrary(tibble)
#sfExport('rot')
#linesswaths <- sfLapply(1:nrow(tr),make_swaths,mod=modquads,emod=eyemod,rocimod=ROCImod,pocimod=POCImod,tracks=tr,minpresss=minpress)
#sfStop()

##  The function returns results in two forms, linestrings and polygons
##  separate into their respective parts
swaths <- do.call(rbind,lapply(linesswaths,'[[',1))
linestrings <-  do.call(rbind,lapply(linesswaths,'[[',2))


###  visualize the result
bind_rows(swaths %>% mutate(geomtype="polys"),linestrings %>%mutate(geomtype="lines"))%>%
arrange(kts)%>%
ggplot()+
geom_sf(data=rnaturalearth::ne_countries(country=c("united states of america","mexico","canada"),scale="medium"),fill=NA)+
geom_sf(aes(col=kts,fill=kts))+
scale_color_gradientn(
    colours = c("black","lightgreen","#44AA99","#DDCC77","#CC6677","#882255"),
    values = scales::rescale(c(64, 83, 96, 113,137)), # breakpoints in data space
    limits = c(0, 150))+
scale_fill_gradientn(
    colours = c("black","lightgreen","#44AA99","#DDCC77","#CC6677","#882255"),
    values = scales::rescale(c(64, 83, 96, 113,137)), # breakpoints in data space
limits = c(0, 150))+
facet_wrap(~geomtype)+
ggtitle("Hurricane Camille wind extents")+
theme_bw()+
xlim(100,80)+ylim(20,45)

```
The maps show the lines and polygons versions of the resulting wind extents. Notice how the extents are quadrant specific. Again, none of these wind extents were provided in the original tabular data, they have been predicted from the models constructed above, mostly based on the provided minimum pressure at each time step. The lower graph singles out one of the time steps around when the storm made landfall in Louisiana on August 17th, 1969. 

![](README_files/figure-gfm/Camille-linesandpolys.png)<!-- -->
![](README_files/figure-gfm/Camille-linesandpolys-landfall.png)<!-- -->
    
These datasets are now similar to what the National Hurricane Center produces for most modern storms, but not for any storms before 2009. Thus, to provide a consistent set of spatial wind extents for any storm before around 2010, we need to use the tabular information provided by IBTrACS and convert it to spatial information as demonstrated above. 

## Interpolate wind and pressure fields

The above steps only provides us with the extent of different wind speeds in 3-hour increments. We need more fine-grained detailed information on what happens between the different wind extents and between the three hour time steps. To achieve this, we combine linear interpolation with a thin plate spline interpolation and arrive at a much smoother representation of wind fields within the storm.

First read in the same tabular date from IBTrACS and convert the individual center points to lines representing the 3 hour track segments of each storm

``` r
###  make into lines from points
###  read in the tabular data
IBTrACS <- read.csv("ibtracs.ALL.list.v04r01.csv") %>%
  filter(SEASON==1969,NAME=="CAMILLE")%>%
  ##  make sure the timestamp is appropriately formatted
  mutate(ISO_TIME = as.POSIXct(ISO_TIME,tz="GMT"),
  BASIN=if_else(is.na(BASIN),"NA",BASIN)) %>%
  ###  turn into point geometry
  st_as_sf(coords=c("LON","LAT"),crs=4326,remove=FALSE) %>%
  ###  the lines are two sets of points, so we need to set the leading point and the lagging point
  mutate(
    geometry_lead = lead(geometry, default = NULL)
  ) %>%
  # drop the NA row created by lagging
  slice(-n()) %>% 
  ###  now combine the two sets of points
  mutate(line = st_sfc(purrr::map2(
      .x = geometry, 
      .y = geometry_lead, 
      .f = ~{st_union(c(.x, .y)) %>% st_cast("LINESTRING")}
    ))) %>%
  st_drop_geometry()%>%
  select(-geometry_lead) %>%
  rename(geometry="line") %>%
  st_set_geometry("geometry") %>%
  st_set_crs(st_crs(4326)) %>% arrange(ISO_TIME)

##  Now get the lines produced earlier and representing the wind extents 
lne <- linestrings %>%
  ## make sure the dates coincide
  filter(date>=min(IBTrACS$ISO_TIME,na.rm=T),date<=max(IBTrACS$ISO_TIME)) %>%
  arrange(date,kts) %>% 
  group_by(date) %>%
  mutate(l = seq(1,n())) %>%
  ungroup() %>%
  ##  get the time difference between the start and end of the line segment
  mutate(dt=as.numeric(difftime(as.POSIXct(date),as.POSIXct(diag(sapply(l, function(n) dplyr::lag(date, n=n))))))) %>%
  st_transform("+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
  
###  interpolate points along a given track segment
###  track segments are generally 3 hours in duration
###  this is too coarse for the raster interpolation
###  this function provides new points along the line in 3 minute intervals
###  it assumes the storm is moving at the same speed during the 3 hour segment
trck_points <- track_interp(IBTrACS) %>%
    st_transform(st_crs(lne))
###  get the date stamps of the original track
###  we will iterate over these
dts = unique(lne$date)[1:(length(unique(lne$date))-1)]

###  Now the 'meat' of the process is to iterate over each line segment, create end-point rasters of wind speed for the beginning and end of the segment 
###  using Thin Plate Spline to interpolate between the lines we produced earlier. 
###  Then linearly interpolate between those end points to get the 3 min interval windspeeds
###  The TPS interpolation is resource hungry and may be limited on insufficent machines. There are alternatives but those arent discussed here.

###  the non-parallel version, much slower, is
fields <- sapply(dts,FUN=fields_wrapper,dates=dts,line=lne,track_points=trck_points,oldw=getOption("warn"),track=IBTrACS,parll=F)

### the parallel version
#library(snowfall)
#sfInit(parallel=TRUE, cpus=10)#, slaveOutfile=logtmp)
#sfLibrary(sf)
#sfLibrary(terra)
#sfLibrary(dplyr)
#sfLibrary(fields)
#sfLibrary("snowfall", character.only=TRUE)
#sfExport('move_and_interp_snowfall')
#sfExport('line_int_snowfall')
#sfExport('comparewinds')
#sfExport('get_dir')
fields <- sfLapply(dts,fun=fields_wrapper,dates=dts,line=lne,track_points=trck_points,oldw=getOption("warn"),track=IBTrACS,parll=T)

```

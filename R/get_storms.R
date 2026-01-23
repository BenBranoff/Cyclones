#' @importFrom dplyr filter mutate any_of select pull slice group_split syms coalesce group_keys
get_storms <- function(source="ncei",id=NULL,name=NULL,season=NULL,basin=NULL,ib_filt=NULL){
  tmf <- tempfile()
  if (is.data.frame(source)){
    dat=source
  }else if (is.character(source)){
    if(source=="ncei") {
      cat("Downloading IBTrACS from: https://www.ncei.noaa.gov/products/international-best-track-archive")
      if (is.null(ib_filt)||isTRUE(as.numeric(season)>=(as.numeric(format(Sys.time(),"%Y"))-3))) ib_filt = "last3years" else if (is.null(ib_filt)||isTRUE(as.numeric(season)>1980)) ib_filt= "since1980"
      ## get appropriate url
      if (!is.null(ib_filt)){
        if (!ib_filt %in% c("ACTIVE","ALL","EP","NA","NI","SA","SI","SP","WP","last3years","since1980")) stop("ib_filt value not found in https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/") else
          URL <- paste0("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.",ib_filt,".list.v04r01.csv")
      }else{
        if (!is.null(basin)){
          if (!basin %in% c("EP","NA","NI","SA","SI","SP","WP")) stop("Provided basin value not found in IBTrACS") else
            URL =  paste0("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.",basin,".list.v04r01.csv")
        }else{
          URL <- paste0("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.ALL.list.v04r01.csv")
        }
      }
      ###  download with url
      tryCatch({
        # Attempt to download the file
        download.file(URL,paste0(tmf,".csv"))
        dat = read.csv(paste0(tmf,".csv"))
        message("Download successful!")
      },
      error = function(e) {
        # Check if the error message indicates a timeout
        if (grepl("Timeout", e$message, ignore.case = TRUE)) {
          warning("Download timeout reached. Please check your internet connection or increase the R timeout option. Alernatively, download and source data directly from browser.")
        } else {
          # Handle other potential errors (e.g., file not found, permission issues)
          warning(paste("Download failed with error:", e$message))
        }
        # Return NULL or an indicator of failure
        invisible(NULL)
      },
      warning = function(w) {
        if (grepl("Timeout|!= reported length", w$message, ignore.case = TRUE)) {
          stop("Download timeout reached. Please check your internet connection or increase the R timeout option. Alternatively, download and source data directly from browser.")
        }else{
          # Handle other warnings, such as "downloaded length != reported length"
          message(paste("A warning occurred during download:", w$message))
        }
      })
      ##  HURDAT
    }else if (source=="hurdat"){
      if (is.null(basin))  stop("Basin not specified. HURDAT available for North Atlantic 'NA' and 'EP' basins.")
      else
        if (basin=="NA"){
          cat("Downloading North Atlantic HURDAT2 from: https://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html")
          url = "https://www.aoml.noaa.gov/hrd/hurdat/hurdat2.html"
        }else if(basin=="NP"){
          cat("Downloading Northeast HURDAT2 from: https://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html")
          url = "https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-nepac.html"
        }else{
          stop("Basin not available in HURDAT. HURDAT available for North Atlantic 'NA' and 'EP' basins.")
        }
      tryCatch({
        # Attempt to download the file
        download.file(url,paste0(tmf,".txt"))
        hurdat <- readLines(paste0(tmf,".txt"))
        message("Download successful!")
      },
      error = function(e) {
        # Check if the error message indicates a timeout
        if (grepl("Timeout", e$message, ignore.case = TRUE)) {
          warning("Download timeout reached. Please check your internet connection or increase the R timeout option. Alernatively, download and source data directly from browser.")
        } else {
          # Handle other potential errors (e.g., file not found, permission issues)
          warning(paste("Download failed with error:", e$message))
        }
        # Return NULL or an indicator of failure
        invisible(NULL)
      },
      warning = function(w) {
        if (grepl("Timeout", w$message, ignore.case = TRUE)) {
          stop("Download timeout reached. Please check your internet connection or increase the R timeout option. Alernatively, download and source data directly from browser.")
        }else{
          # Handle other warnings, such as "downloaded length != reported length"
          message(paste("A warning occurred during download:", w$message))
        }
      })
      start <- grep("<body>",head(hurdat)) +2
      hurdat <- hurdat[start:length(hurdat)]
      ls <- lengths(strsplit(hurdat,","))
      hdlines <- which(ls==ls[1])
      idlines <- hurdat[hdlines]
      stormids <- sapply(strsplit(idlines,","),"[[",1)
      stormlengths <- as.numeric(sapply(strsplit(idlines,","),"[[",3))
      RLES <- if (is.null(id)) cbind(hdlines, stormlengths) else cbind(hdlines[grep(paste(id,collapse="|"),stormids)],stormlengths[grep(paste(id,collapse="|"),stormids)])
      RLES <- split(RLES,seq(nrow(RLES)))
      storms <- lapply(RLES,function(x,y){
        read.table(text=y[(x[1]+1):(x[1]+x[2])],sep=",",header=FALSE,
                   col.names=c("date","time","recID","Status","LAT","LON","USA_MSW","USA_PRES",
                               "USA_R34_NE", "USA_R34_SE","USA_R34_SW","USA_R34_NW",
                               "USA_R50_NE", "USA_R50_SE","USA_R50_SW","USA_R50_NW",
                               "USA_R64_NE", "USA_R64_SE","USA_R64_SW","USA_R64_NW","USA_RMW")) |>
          mutate(USA_ATCF_ID=strsplit(y[x[1]],",")[[1]][1],
                 NAME=gsub("^\\s+|\\s+$","",strsplit(y[x[1]],",")[[1]][2]),
                 ISO_TIME=as.POSIXct(paste0(substr(date,1,4),"-",substr(date,5,6),"-",substr(date,7,8)," ",sprintf("%04d",time)),format="%Y-%m-%d %H",tz="UTC"))
      },y=hurdat)
      dat <- do.call(rbind,storms) |>
        mutate(BASIN=basin,SEASON=substr(date,1,4))
      warning("HURDAT data missing ROCI column. This value will be modelled when needed to provide maximum storm extent.")
    } else {
      if (!grepl(".csv",source)) stop("file must be sourced from 'ncei', 'hurdat' or from a .csv")
      dat=read.csv(source)
    }
  }else{
    stop("unknown data source")
  }
  if (any(grepl("nmile",dat|>slice(1)))) dat <- dat |> slice(-1)
  dat=dat |>
    select(any_of(c("SID","USA_ATCF_ID","SEASON","NAME","BASIN","ISO_TIME","LAT","LON","USA_WIND","USA_PRES","USA_SSHS",
                  "USA_R34_NE", "USA_R34_SE","USA_R34_SW","USA_R34_NW",
                  "USA_R50_NE", "USA_R50_SE","USA_R50_SW","USA_R50_NW",
                  "USA_R64_NE", "USA_R64_SE","USA_R64_SW","USA_R64_NW",
           "USA_ROCI","USA_POCI","USA_RMW","REUNION_RMW","BOM_RMW","USA_EYE","STORM_SPEED"))) |>
    mutate(ISO_TIME=as.POSIXct(ISO_TIME,tz="UTC"),
           BASIN=if_else(is.na(BASIN),"NA",BASIN))
  if (!is.null(basin)){
    if (toupper(basin) %in% c("EP","NA","NI","SA","SI","SP","WP"))
      dat <- dat |> filter(BASIN %in% toupper(basin))
    else
      stop("Provided basin value not found in data. Check source and other filters")
  }
  if (!is.null(season)){
    if (as.numeric(season) %in% seq(1851,format(Sys.time(),"%Y")))
      dat <- dat |> filter(as.numeric(season) %in% SEASON)
    else
      stop("Provided season value not found in data. Check source and other filters.")
  }
  if (!is.null(name)){
    names=dat |> pull(NAME)|>unique()
    if (sum(grepl(paste0(toupper(name),collapse="|"),names),na.rm=TRUE)>0){
      if (sum(grepl(paste0(toupper(name),collapse="|"),names),na.rm=TRUE)<length(name)){warning("Not all names found in data. Check other filters.")}
      dat <- dat |> filter(NAME %in% toupper(name))
    }else{ stop("Provided season value not found in data. Check source and other filters.")}
  }
  if (!is.null(id)){
    if ((id %in% unique(dat$SID)|id %in% unique(dat$USA_ATCF_ID)))
      dat <- dat |> filter(SID %in% id|USA_ATCF_ID %in% id)
    else
      stop("Provided id value not found in data. Check source and other filters. HURDAT uses the USA_ATCFID while IBTrACS uses either USA_ATCFID or SID.")
  }
  dat <- dat |>
    mutate(ID=coalesce(!!!syms(intersect(c("SID","USA_ATCF_ID"),names(dat)))),
           ID=paste(NAME,SEASON,BASIN,ID, sep="_")) |>
    tidyr::nest(.by =ID)
  names(dat$data) <- dat$ID
  #ids <- dat |>
  #  group_keys() |>
  #  pull(ID)
  #dat <- dat |> group_split() |>
  #  setNames(ids)
  return(dat)
}


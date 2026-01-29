checktracks <- function(dta){
  options(warn = 1)
  if (paste0(class(dta),collapse="")=="vctrs_list_ofvctrs_vctrlist"&&length(dta)>1){
    warning("More than one storm in data. Only first will be used. Use lapply or a parallel equivalent to repeat for multiple storms.")
    return(dta[[1]]) 
  } else if(paste0(class(dta),collapse="")=="tbl_dftbldata.frame") { return(dta)}else{
    stop("Unknown data input. Use get_storms to retrieve and/or format data for further processing.")
  }
}
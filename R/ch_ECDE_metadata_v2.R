#'
#'  new version of ch_ECDE_metadata
#'
#'


ch_get_ECDE_metadata_v2 <- function (filename, writefile = NULL)
{
  if (filename == "" | is.null(filename)) {
    stop("ECDE file not specified")
  }
  if (!file.exists(filename)) {
    stop("ECDE file not found")
  }

  meta <- read.table(filename, skip = 103, sep = " ", na.strings = -999)

  names(meta) <- c("Station", "Fav", "StationName", "HydStatus",
                   "Prov", "Latitude", "Longitude",
                   "DrainageArea", "Eff_DrainageArea", "Years",
                   "From", "To", "Reg.", "Flow", "Level", "Sed", "OperSched",
                   "RealTime", "RHBN", "Region", "Datum", "Operator")
  meta <- meta[, c(1, 3:22)]
  if (!is.null(writefile))
    write.csv(meta, writefile, row.names = FALSE)
  return(meta)
}

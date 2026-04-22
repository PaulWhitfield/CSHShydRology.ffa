library(CSHShydRology)
require(pracma)
library(ppclust)
library(odetector)
library(fcvalid)
require(rrcov)
require(POT)

source("C:/R-packages/ffaScreening/R/S&R_functions testing version.R")
source("C:/R-packages/ffaScreening/R/ffa_screen_plot  reorganized.R")
source("C:/R-packages/ffaScreening/R/ch_ECDE_metadata_v2.R")
source("C:/R-packages/ffaScreening/R/ch_high_Grubbs_test.R")
source("C:/R-packages/ffaScreening/R/mretlev.uvplot for POT.R")
source("C:/R-packages/ffaScreening/R/ch_circular_colors.R")

options(warn = 1)
setwd("c:/p_flood_of_record")

meta <- ch_get_ECDE_metadata_v2("FavHydatStations_dhb.tb0",
                                writefile = "dhb flood.csv")

setwd("c:/p_flood_of_record/dhb_list")

pdf("ffa_screen plot POT with timing and polar versions 2.pdf")

print("Testing Pareto axis corrected from years to events")

files <- list.files(pattern = "_ts.csv")


for( kk in c(37, 1, 4,38)){

  mdata <- ch_read_ECDE_flows(files[kk])

### use minimum amax as threshold
  n0peaks <- ch_sh_get_amax(mdata)
  n0peaks <- n0peaks[n0peaks$days >= 365,]  # remove partial years
  thresh <- min(n0peaks$amax)

### get peaks over threshold
  npeaks <- ch_get_peaks(mdata,threshold = thresh)
  mpeaks <- npeaks$POTevents

### get station metadata
  stn <- substr(files[kk],1, 7)

  station <- meta[meta$Station == stn,]
  print(paste(kk, stn, "#events=", length(mpeaks[,1])))
  mtitle <- paste(station$Station,station$StationName,"-",station$Prov)
  par(oma = c(1,1,1,1))

  mplot0 <- ch_ffa_screen_plot(mpeaks, mtitle = mtitle,
              mcol = c("black","gray40","black","purple"),stn = stn)


}

 graphics.off()

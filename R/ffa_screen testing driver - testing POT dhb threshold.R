library(CSHShydRology)
require(pracma)
library(ppclust)
library(odetector)
library(fcvalid)
require(rrcov)
require(POT)
require(DescTools)
require(tscount)

source("C:/R-packages/ffaScreening/R/SR_hstat.R")
#source("C:/R-packages/ffaScreening/R/ffa_screen_plot  reorganized edits.R")
source("C:/R-packages/ffaScreening/R/ch_ECDE_metadata_v2.R")
source("C:/R-packages/ffaScreening/R/ch_high_Grubbs_test.R")
source("C:/R-packages/ffaScreening/R/ch_plot_POT.R")
source("C:/R-packages/ffaScreening/R/ch_circular_colors.R")
source("C:/R-packages/ffaScreening/R/rev ch_get_peaks.R")

graphics.off()

options(warn = 1)
setwd("c:/p_flood_of_record")

meta <- ch_get_ECDE_metadata_v2("FavHydatStations_dhb.tb0",
                                writefile = "dhb flood.csv")

setwd("c:/p_flood_of_record/dhb_list")

pdf("ffa_screen plot POT with timing and polar versions 8b.pdf")

files <- list.files(pattern = "_ts.csv")

thresh0 <- array(0, dim = 40)
  thresh0[1] <- 1000
  thresh0[4] <- 40
  thresh0[37] <- 155
  thresh0[38] <- 445


for( kk in c(37, 1, 4,38)){


  mdata <- ch_read_ECDE_flows(files[kk])

### use minimum amax as threshold


### get peaks over threshold
  npeaks <- ch_get_peaks(mdata,threshold = thresh0[kk])

  mpeaks <- npeaks$POTevents

### get station metadata
  stn <- substr(files[kk],1, 7)

  station <- meta[meta$Station == stn,]
  cat("\n")
  print(paste(kk, stn, "#events=", length(mpeaks[,1])))
  mtitle <- paste(station$Station,station$StationName,"-",station$Prov)
  par(oma = c(1,1,1,1))

  print(paste("POT threshold =", thresh0[kk]))

  mplot0 <- ch_ffa_screen_plot(mpeaks,
                               mtitle = mtitle, thresh = thresh0[kk],stn = stn)


}

 graphics.off()

#'  update ch_get_peaks
#'
#'
ch_get_peaks <- function (dataframe, threshold)
{
  maxflow <- max(dataframe$Flow)
  if (maxflow < threshold) {
    message(paste("Threshold of", threshold, "\n                is greater than maximum observed flow",
                  maxflow))
    return()
  }
  data <- dataframe$Flow
  Date <- dataframe$Date
  inSYM <- dataframe$SYM  ##added
  event <- array(0, dim = length(data))
  event_num <- array(0, dim = length(data))
  flow <- array(dim = 7)
  st_date <- array(NA, dim = 3)
  max_date <- array(NA, dim = 3)
  class(st_date) <- "Date"
  class(max_date) <- "Date"
  case <- list(dim = 3)
  SYM <- array(NA, dim = 3)  ### added
  for (i in 1:length(data)) {
    if (is.na(data[i]))
      next
    if (data[i] > threshold)
      event[i] <- 1
  }
  max <- array(NA, dim = 1)
  volume <- array(NA, dim = 1)
  duration <- array(NA, dim = 1)
  index = 1
  flag = 0
  for (k in 1:length(data)) {
    if (event[k] == 1 && flag == 0) {
      st_date[index] <- Date[k]
      max[index] = data[k]
      max_date[index] <- Date[k]
      SYM[index] <- inSYM[k] ## added
      volume[index] <- 0
      duration[index] <- 0
      flag <- 1
    }
    if (event[k] == 1 && flag == 1) {
      if (data[k] > max[index]) {
        max[index] <- data[k]
        max_date[index] <- Date[k]
        SYM[index] <- inSYM[k] ## added
      }
      volume[index] <- volume[index] + data[k]
      duration[index] <- duration[index] + 1
    }
    if (event[k] == 0 && flag == 1) {
      index <- index + 1
      flag <- 0
    }
  }
  st_date <- as.Date(st_date, format = "%Y-%m-%d")
  volume <- volume * 24 * 60 * 60 * 1e-09
  max_date <- as.Date(max_date, format = "%Y-%m-%d")


  POT_events <- data.frame(st_date, max_date, max, volume,
                           duration, SYM)  ## note change
  flag = 0
  index = 1
  flow <- array(dim = 9)
  for (k in 1:length(data)) {
    if (event[k] == 1 && flag == 0) {
      st_date[index] <- as.character(Date[k])
      if (k == 1) {
        flow[1:4] <- NA
      }
      if (k == 2) {
        flow[1:3] <- NA
        flow[4] <- data[k - 1]
      }
      if (k == 3) {
        flow[1:2] <- NA
        flow[3] <- data[k - 2]
        flow[4] <- data[k - 1]
      }
      if (k == 4) {
        flow[1] <- NA
        flow[2] <- data[k - 3]
        flow[3] <- data[k - 2]
        flow[4] <- data[k - 1]
      }
      if (k <= 4) {
        event_num[1:4] <- index
      }
      if (k > 5) {
        event_num[k - 1] <- index
        flow[4] <- data[k - 1]
        event_num[k - 2] <- index
        flow[3] <- data[k - 2]
        event_num[k - 3] <- index
        flow[2] <- data[k - 3]
        event_num[k - 3] <- index
        flow[1] <- data[k - 4]
      }
      ii <- k
      ii <- ii - 5
      event[k - 1] <- 1
      event_num[k] <- index
      flag <- 1
    }
    if (event[k] == 1 && flag == 1) {
      event_num[k] <- index
      flow[k - ii] <- data[k]
    }
    if (event[k] == 0 && flag == 1) {
      event_num[k] <- index
      flow[k - ii] <- data[k]
      if ((k - length(data)) >= 1)
        event_num[k + 1] <- index
      flow[k + 1 - ii] <- data[k + 1]
      if ((k - length(data)) >= 2)
        event_num[k + 2] <- index
      flow[k + 2 - ii] <- data[k + 2]
      if ((k - length(data)) >= 3)
        event_num[k + 3] <- index
      flow[k + 3 - ii] <- data[k + 3]
      case[index] <- list(flow[!is.na(flow)])
      rm(flow)
      flow <- array(dim = 9)
      event[k - 1] <- 1
      index <- index + 1
      flag <- 0
    }
  }
  ncases <- length(POT_events$st_date)
  events <- list(POT_events, ncases, case)
  names(events) <- c("POTevents", "ncases", "case")
  return(events)
}



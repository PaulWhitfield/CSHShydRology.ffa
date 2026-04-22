#'  revised ch_sh_get_amax
#'
#'

ch_get_amax <- function (df)
{

  data <- df$Flow
  Date <- df$Date
  df$SYM[df$SYM == " "] <- "" #added
  df$SYM[df$SYM == "  "] <- "" #added

  year <- format(Date, "%Y")
  Year <- as.numeric(unique(year))
  maxdate <- array(NA, dim = length(Year))
  doy <- array(NA, dim = length(Year))
  days <- array(NA, dim = length(Year))
  SYM <- array(NA, dim = length(Year))
  class(maxdate) <- "Date"
  year <- as.factor(year)
  amax <- as.numeric(tapply(data, year, max))
  dataframe <- data.frame(df, year)
  for (k in 1:length(Year)) {
    ndata <- dataframe[dataframe$year == Year[k], ]
    days[k] <- length(ndata$Flow)
    ndata <- ndata[ndata$Flow == amax[k], ]
    SYM[k] <- ndata[1, 5]  ### added to get SYM
    maxdate[k] <- ndata[1, 3]
    maxdate_a <- timeDate::as.timeDate(maxdate[k])
    doy[k] <- timeDate::dayOfYear(maxdate_a)
  }
  result <- data.frame(Year, amax, maxdate, doy, days, SYM)
  return(result)
}

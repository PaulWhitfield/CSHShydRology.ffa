#'  ch_ffa_screen_plot
#'
#' generate a flood frequency plot with symbols indicating whether an observation is a
#' high or low outlier, and colour symbols indicating month of the year in which the flood
#' occurred.
#'
#' @param df dataframe with date of event (maxdate) and annual maximum (amax)
#' @param mtitle title for plot
#' @param stn stationID
#' @param trans transform series peak. options are "none", "log", and "boxcox"
#' @param polar produces a polar plot if polar = TRUE default.
#' @param thresh POT threshold default is NULL
#' @param kthresh kappa threshold default is 5
#' @param mcol array of three colours, default is c("orange", "gray70", "red") use to emphasize
#' low, normal, and high values
#' @param n default is 12, number of colours for number of months
#' @param m smallest angle in radians parameter for generating circular colours default = 0
#' @param M largest angle in radians parameter for generating circular colours default = 2 (*pi)
#' @param offset the zero in radians, default is 0.
#' @param polar if TRUE [default] a polar plot of the annual maxima and outliers is produced.
#'
#' @return a list containing
#' \describe{
#'   \item (station) (stationID)
#'   \item (nyears) (number of years of data)
#'   \item (maxQ) (maximum flow - Flood of Record)
#'   \item (mindex) (year with maximum flow)
#'   \item (maxDate) (Date of Flood of Record)
#'   \item (maxdoy) (day of year with FoR)
#'   \item (WalfWolfowitz) (result of WW test against trend)
#'   \item (MannKendall) (result of Mann Kendall trend test)
#'   \item (Pettitt) (result of Pettitt test for single change point)
#'   \item (POTmodel) (model of POT)
#'   \item (out_hgrubb) (number of high outliers, indexes in df)
#'   \item (out_lgrubb) (number of low outliers, indexes in df)
#'   \item (out_timing) (number of timing outliers, indexes in df)
#'   \item (fz_out) (indexes of fuzzy clusters outliers, indexes in df)
#'   \item (df) (a dataframe with results for each year)
#'  }
#' Also produces an annotated three plots: an annotated time series plot,
#' an return period plot of flood events and an optional polar plot
#'
#' @references
#'
#' Cohn, T. A., J. F. England, C. E. Berenbrock, R. R. Mason, J. R.
#' Stedinger and J. R. Lamontagne (2013). "A generalized Grubbs‐Beck test statistic
#' for detecting multiple potentially influential low outliers in flood series."
#' Water Resources Research 49(8): 5047-5058 10.1002/wrcr.20392: 10.1002/wrcr.20392.
#'
#' Sau, M. F. and D. Rodriguez (2018). "Minimum distance method for directional data and outlier detection." Advances in Data Analysis and Classification 12: 587-603 DOI: 10.1007/s11634-017-0287-9.
#'
#' Whitfield, P. H. and D. H. Burn (2025*). "Extreme and Rogue Floods in North America."
#' Journal of Hydrology. in review.
#'
#'
#' @importFrom pracma rad2deg
#' @import circular
#' @import plotrix
#' @import CSHShydRology
#' @importFrom MGBT MGBT
#' @importFrom outliers grubbs.test
#' @importFrom DescTools IsLeapYear
#' @importFrom Hmisc subplot
#'
#' @export
#'@examples \donttest{
#'# Not tested automatically as can be very slow to execute
#' data(CAN05AA008)
#' amax <- ch_get_amax(CAN05AA008)
#' ch_ffa_screen_plot(amax)}



ch_ffa_screen_plot <- function (df, mtitle = NULL, stn,
                                polar = TRUE,
                                kthresh = 5,
                                thresh = NULL,
                                n = 12, m = 0, M = 2, offset = 0,
                                trans = "log",
                      mcol = c( "orange","gray60", "red", "darkred"))
  {


  ##### determine type of data
  if(names(df)[1] == "Year"){
    print("data is amax")
    POT = FALSE
    AMAX = TRUE
    fitted = NULL  ## needs to exist for POT result to be available
    df$c0 <- rep(2, length(df$amax))
    df$c0[df$days <= 364] <- 1

  }

  if(names(df)[1] == "st_date"){
    #### this is POT data
    print(paste("data is POT - threshold = ", thresh))
    POT = TRUE
    AMAX = FALSE
    df0 <- df
    maxdate <- df$max_date
    amax <- df$max
    SYM <- df$SYM
    maxdate_a <- timeDate::as.timeDate(df$max_date)
    doy <- timeDate::dayOfYear(maxdate_a)
    Year <- as.numeric(format(maxdate, "%Y"))
    days <- array(365, length(Year))

   c0 <- rep(2, length(df$max))

    for (jj in 1 : length(Year)){
    if(DescTools::IsLeapYear(Year[jj])) days[jj] <- 366
    }

  df <- data.frame(Year,amax,maxdate,doy,days, SYM, c0)
  }

  mnths <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

  mon <- as.numeric(format(df$maxdate, "%m"))
  df$mon <- mon

  m3s <- expression(paste("m"^3,"/s"))



  mcode <-array(2, dim = length(df$amax))

  nyears <- length(df$amax)
  maxQ <- max(df$amax, na.rm = TRUE)

  mindex <- which.max(df$amax)

  maxDate <- df$maxdate[mindex]

  maxdoy <- df$doy[mindex]

  mmax <- sort(df$amax, decreasing = TRUE)

  mmax<- mmax[mmax < maxQ]

fpch <- c(1,19)

##########################################  Part 1 time series plot
  par(mar = c(3,5,3,5))

  med_amax <- median (df$amax)
 xlims <- c(min(df$Year), max(df$Year))
 ylims <- c(min(df$amax), max(df$amax))

  plot(df$Year[df$amax == med_amax], df$amax[df$amax == med_amax], main = mtitle,
       xlim = xlims, ylim = ylims,
       pch = fpch[df$c0], col = "gray30", cex = 0.8,
       ylab = expression( "Peak Flow ("*m^3*s^{-1}*")" ) ,
       xlab = "", las = 1)

  abline(h = med_amax, col = "grey50", lty = 3)

  s <- sign(df$amax - med_amax)
  for (i in 1:(length(df$amax) - 1)) {
    if (s[i] * s[i + 1] < 0) {
      abline(v = df$Year[i] + 0.5, lty = 2, col = "grey70")
    }
  }

  points(df$Year[df$amax > med_amax], df$amax[df$amax > med_amax], pch = fpch[df$c0], col = "blue", cex =0.8)
  points(df$Year[df$amax < med_amax], df$amax[df$amax < med_amax], pch = fpch[df$c0], col = "red" ,cex = 0.8)

 ##################  add symbols to points flagged by WSC
    df$SYM[df$SYM == "A"] <- 1
    df$SYM[df$SYM == "B"] <- 2
    df$SYM[df$SYM == "D"] <- 3
    df$SYM[df$SYM == "E"] <- 4
    df$SYM <- as.numeric(df$SYM)
    sym_col <- c("green3", "cyan3", "orange", "red")
    points (df$Year, df$amax, pch = 22, cex = 1.3, font = 2,col = sym_col[df$SYM])
 ############################### end flagging
 par(xpd = TRUE)

 legend("right", legend = c("pArtial day", "Backwater", "Dry", "Estimate"), pch = 22,
        col = sym_col, cex = 0.55 , title = "SYM", inset = c(-0.17,0.0))

if(AMAX)
  legend("bottomright", legend = c("<365", "365"), pch = fpch,
        col = "black", cex = 0.55 , title = "Days", inset = c(-0.12,0.0))

 par(xpd = FALSE)




  ############################################# Mann Kendall
  mk <- Kendall::MannKendall(df$amax)

  trend_line <- predict(loess(df$amax ~ df$Year))
  if(mk$sl > 0.1)
    lines(df$Year, trend_line, lty = 2, col =  "gray30", lwd = 1.75)
  if(mk$sl <= 0.1)
    lines(df$Year, trend_line, lty = 2, col =  "darkred", lwd = 1.75)
  if(mk$sl <= 0.05)
    lines(df$Year, trend_line, lty = 1, col =  "red", lwd = 2.5)

  ###################################### Wald-Wolfowitz against trend
  ## left.sided tests against trend
  wwo <- randtests::runs.test(df$amax, pvalue = "normal",
                              alternative = "left.sided",
                              plot = FALSE)

  ########################################## missing years
  pet <- trend::pettitt.test(df$amax)
  if(pet$p.value > 0.10) abline (v = df$Year[pet$estimate], col = "gray30", lty = 3,lwd = 1.4)
  if(pet$p.value <= 0.10) abline (v = df$Year[pet$estimate], col = "darkred", lty = 2, lwd = 1.75)
  if(pet$p.value <= 0.05) abline (v = df$Year[pet$estimate], col = "red", lty = 1, lwd = 2.5)

 ########################################## missing years
  myear <- c(min(df$Year): max(df$Year))
  missyear <- myear[!myear %in% df$Year]


  for (jj in 1: length(missyear)){
    abline( v = missyear, col = "pink2", lwd = 1.7)
  }


  ################################# add summary information
if(length(missyear) == 0)
mtext(paste ("#Years without obs =",length(missyear)),
          side = 3, line = -1, col = "black", cex = 0.75, adj = 0.03, font = 1)
if(length(missyear) >= 1)
    mtext(paste ("#Years without obs =",length(missyear)),
          side = 3, line = -1, col = "red", cex = 0.75, adj = 0.03, font = 2)

if (pet$p.value > 0.05)
  mtext(paste ("Pettitt changepoint",df$Year[pet$estimate],"p value", round(pet$p.value, 3)),
        side = 3, line = -1, col = "black", cex = 0.75, adj = 0.9)
if (pet$p.value <= 0.05)
    mtext(paste ("Pettitt changepoint",df$Year[pet$estimate],"p value", round(pet$p.value, 3)),
          side = 3, line = -1, col = "red", font = 2, cex = 0.75, adj = 0.9)


if(wwo$p.value > 0.05)
  mtext(paste("WW runs test against trend p value",round(wwo$p.value,3)), side = 1,
        line = -1, cex = 0.7, adj = 0.1)
if(wwo$p.value <= 0.05)
    mtext(paste("WW runs test against trend p value",round(wwo$p.value,3)), side = 1,
          line = -1, cex = 0.7, adj = 0.1, font = 2, col = "red")

if(mk$sl > 0.05)
  mtext(paste("Mann-Tendall tau", round(mk$tau, 3), "p value", round(mk$sl,3)),
        side = 1, line = -1, cex = 0.7, adj = 1.0)
if(mk$sl <= 0.05)
  mtext(paste("Mann-Tendall tau", round(mk$tau, 3), "p value", round(mk$sl,3)),
        font = 2, col = "red",
          side = 1, line = -1, cex = 0.7, adj = 1.0)

  if(POT) mtext(paste("POT events threshold=", thresh),
                side = 3, line = 0.02, adj =0.05, cex = 0.65)


##################################################### time series plot end




############### set variables for case of no fuzzy clustering
   cl_num <- array( 1, nyears)
   mcases <- NA
##########################################  Part 2 seek outliers
######################### get codes Grubbs

############### high Grubbs


#### revised after HSJ review



if(trans == "none"){
 print("using untransformed values")
 gt0 <- outliers::grubbs.test(df$amax, type = 10, opposite = FALSE)
 gtest <- ch_high_Grubbs_test(df$amax)
if(gt0$p.value <= 0.05 && substr(gt0$alternative,1,5) == "lowes") gtest$tout[gtest$tout == 1] <- 0
}

if(trans == "log"){
  print("using log transformed values")
 gt0 <- outliers::grubbs.test(log(df$amax), type = 10, opposite = FALSE)
 gtest <- ch_high_Grubbs_test(log(df$amax))
if(gt0$p.value <= 0.05 && substr(gt0$alternative,1,5) == "lowes") gtest$tout[gtest$tout == 1] <- 0
}

if(trans == "boxcox"){
  print("using boxcox transformed values")
  gdata <- bcTransform(df$amax, plot = FALSE, verbose = FALSE)$tf.data
  gt0 <- outliers::grubbs.test(gdata, type = 10, opposite = FALSE)
  gtest <- ch_high_Grubbs_test(gdata)
  if(gt0$p.value <= 0.05 && substr(gt0$alternative,1,5) == "lowes") gtest$tout[gtest$tout == 1] <- 0
}



for( ll in 1: length(mcode)) {
  if(gtest$tout[ll] == 1) mcode[ll] <- 3
}

df$mcode <- mcode

  ############################## low Grubbs
  mg <- MGBT::MGBT(df$amax)

  gindex <- which(df$amax < mg$LOThresh)
  l_grubb <- rep(1,length(df$amax))
  l_grubb[l_grubb < mg$LOThresh] <-2

  for( ll in 1: length(gindex)) {
    mcode[gindex[ll]] <- 1
  }


  ############################## low Grubb end

  c_days <- df$doy/df$days *360
  rad_days <- deg2rad(c_days)

###############  implement Sau & Rodrigues 2018

rad_days <- circular::circular(rad_days, units = "radians",
                                    zero = 3*pi/2, rotation = "clock")



SR_test <- SR_hstat(rad_days)

  isout <- SR_test$isout

 isout[isout == 0] <- 2
 df$isout <- isout


porcout <- SR_test$porcout


#################################### seeking outliers end

############################################# common plotting info

qcol = c("darkblue",  "blue", "red", "orange", "red",
         "darkblue", "darkgreen", "magenta" , "purple", "cyan", "black")
qpch = c(1, 19,  24, 22, 15, 11, 12, 19)
qcex = c(0.7, 1.1, 1.1, 1.1, 1.1, 1.0, 1.1, 1.1, 1.1, 0.6)

sr_pch <- c( NA, 8)
mpch <- c(22, 21, 24, 8 )
mcex <- c(1.20, 0.9, 1.20)

ccol <- ch_circular_colors(n=12, m = m, M = M*pi)
 par(xpd = FALSE)


if(AMAX){
########################################### plotting Gumbel  for amax
# Generate plotting positions

  Q <- df$amax
  n = length(Q)
  r = n + 1 - rank(Q)  # highest Q has rank r = 1
  T = (n + 1)/r

  mat0 <- pretty(0:maxQ)




  # Set up x axis tick positions and labels
  Ttick = c(1.001,1.01,1.1,1.5,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45,50,60,70,80,90,100,200, 300)
  xtlab = c(1.001,NA,NA,NA,2,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,15,NA,NA,NA,NA,20,NA,30,NA,NA,NA,50,NA,NA,NA,NA,100, 200, 300)
  y = -log(-log(1 - 1/T))
  ytick = -log(-log(1 - 1/Ttick))
  xmin = min(min(y),min(ytick))
  xmax = max(ytick)


  # Fit a line by method of moments, along with 95% confidence intervals
  KTtick = -(sqrt(6)/pi)*(0.5772 + log(log(Ttick/(Ttick-1))))
  QTtick = mean(Q) + KTtick*sd(Q)
  nQ = length(Q)
  se = (sd(Q)*sqrt((1+1.14*KTtick + 1.1*KTtick^2)))/sqrt(nQ)
  LB = QTtick - qt(0.975, nQ - 1)*se
  UB = QTtick + qt(0.975, nQ - 1)*se
  ML <- (LB+UB)/2
  max = max(UB)
  Qmax = max(QTtick)

  par(oma = c(1,1,1,1))
  par(mar = c(4,5,3,1))

    # Plot peak flow series with Gumbel axis
    plot(jitter(y,8), jitter(Q,8),
         ylab = expression( "Peak Flow ("*m^3*s^{-1}*")" ) ,
         xaxt = "n", xlab = "Return Period (yr) [Gumbel]",
         col.lab = "black",
         las  = 1,
         tcl = -0.35,
         ylim = c(0, max(Q)),
         xlim = c(xmin, xmax),
         cex = mcex[mcode],
         pch = mpch[mcode],
         col = mcol[mcode],
         bg = ccol[mon],
         main = "" )

# add timing outliers


    points(y,Q, pch = sr_pch[isout], cex = mcex[mcode]*1.1,
           font = 2, col = mcol[4])


    if (nchar(mtitle) <= 50) title(main = mtitle)
    if (nchar(mtitle) > 50) title(main = mtitle, cex.main = 0.95)

    axis(1, at = ytick, labels = as.character(xtlab), cex = 0.65)


    # Add fitted line and confidence limits
    lines(ytick, ML, col = "black")
    lines(ytick, LB, col = "black", lty = 2)
    lines(ytick, UB, col = "black", lty = 2)


      mline <- lm(QTtick ~ ytick)
    abline(h = mat0, col="gray50", lty = 3)

    vlines <- c(ytick[1], ytick[5], ytick[8], ytick[13], ytick[23], ytick[29], ytick[34], ytick[35])
    abline(v = vlines , col="gray50", lty = 3)

    text(0.5, 0, paste( length(Q), "events"), pos=2, cex =0.75)

########### add legend with counts
    ltext <- c(paste("low Grubbs n=", mg$klow),
               paste("normal n=",length(df$amax)-mg$klow-sum(gtest$tout)),
               paste("high Grubbs n=",sum(gtest$tout)),
               paste("Timing outliers n=",length(isout[isout==2])))

################## add subplot as legend
    mnths <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

###############################################################
    legend = "norm"
    if (max(df$amax) <= 350) legend = "flip"

    if (legend == 'norm') {

      legend("topleft", ltext, col = mcol, pch = mpch,
             pt.cex = mcex, cex = 0.6, bg = "white")

      ################## add subplot as legend
      msub <- Hmisc::subplot(hist(df$mon, breaks = c(0:12), col = ccol, xaxt = "n",
                                  bg = "white",
                                  xlab = "", main = "", cex.axis = 0.45,
                                  ylab = "Proportion",
                                  freq = FALSE, las = 1),
                             x = 2.50,  y = 0.0, hadj = 0, vadj = 0,
                             size = c(1.5, 0.75))

      op <- par(no.readonly=TRUE)
      par(msub)
      box()
      axis(1, at = 0.5:11.5, labels = mnths, tick = TRUE, lwd = -1,
           line = -1, cex.axis = 0.6)
      par(op)
    }

    if(legend == "flip"){

      legend("bottomright", legend = ltext, col = mcol,
             pch = mpch, pt.cex = mcex, cex = 0.6, bg = "white")
      par(bg = "white")

      msub <- Hmisc::subplot(hist(df$mon, breaks = c(0:12), col = ccol, xaxt = "n",
                                  bg = "white",
                                  xlab = "", main = "", cex.axis = 0.45,
                                  ylab = "Proportion",
                                  freq = FALSE, las = 1),
                             x = 0.00,  y = maxQ, hadj = 0, vadj = 1,
                             size = c(1.5, 0.75))
      op <- par(no.readonly=TRUE)
      par(msub)
      box()
      axis(1, at = 0.5:11.5, labels = mnths, tick = TRUE, lwd = -1,
           line = -1, cex.axis = 0.6)
      par(op)
    }

 ########
}
############################# end of amax plot

########################################################### POT plot
if(POT){



    fitted <- fitgpd(df[,"amax"], threshold = thresh, est = "mle")
     npy <- length(df$amax)/length(unique(df$Year))

    df1 <- df[order(df$amax),]

    par(oma = c(1,1,1,1))
    par(mar = c(4,5,3,1))

  ############################################ plotting Pareto
  #### ch_plot_POT is mostly from a plot in the package POT
  #### it has been modified to plot different symbols, and different colours
  #### These are mon (month code), mcode (a symbol type 1 2, or 3),
  #### ccol (12 month colour code),  mcol (colours for line around symbol),
  #### isout (an index regarding being normal 1 or a timing outlier 2),
  #### sr_pch (symbols for isout, currently NA and 8)
  #### mcex (main symbol sizes), and  mpch (main symbols 22, 21, 24 for normal and outliers)

      mpareto_test <- ch_plot_POT(fitted,
                           ylab = expression( "Peak Flow ("*m^3*s^{-1}*")"),
                           ci = TRUE, points = TRUE, las = 1,
                           npy = npy,
                           mon = df1$mon,
                           mcode = df1$mcode,
                           ccol = ccol,
                           mcol = mcol,
                           xlimsup = 300.,
                           isout = df1$isout,
                           sr_pch = sr_pch,
                           mcex = mcex,
                           mpch = mpch,
                           main = mtitle)



    mtext(paste("POT threshold = ", thresh, "#py =", round(npy, 3)),side = 1, line = -1, adj =0.1, cex = 0.6)

    #############################################################

    ########### add legend with counts
    ltext <- c(paste("low Grubbs n=", mg$klow),
               paste("normal n=",length(df$amax)-mg$klow-sum(gtest$tout)),
               paste("high Grubbs n=",sum(gtest$tout)),
               paste("Timing outliers n=",length(isout[isout==2])))

    ###############################################################

    legend = "norm"
    if (max(df$amax) <= 350) legend = "flip"

    if (legend == 'norm') {

      legend("topleft", ltext, col = mcol, pch = mpch,
             pt.cex = mcex, cex = 0.65, bg = "white")


      ################## add subplot as legend
      msub <- Hmisc::subplot(hist(df$mon, breaks = c(0:12), col = ccol, xaxt = "n",
                                  bg = "white",
                                  xlab = "", main = "", cex.axis = 0.45,
                                  ylab = "Proportion",
                                  freq = FALSE, las = 1),
                             x = 1.7, y = (1.1 * thresh), hadj = 0, vadj = 0,
                             size = c(1.5, 0.75))

      op <- par(no.readonly=TRUE)
      par(msub)
      box()
      axis(1, at = 0.5:11.5, labels = mnths, tick = TRUE, lwd = -1,
           line = -1, cex.axis = 0.6)
      par(op)
    }

    if(legend == "flip"){

      legend("bottomright", legend = ltext, col = mcol,
             pch = mpch, pt.cex = mcex, cex = 0.6, bg = "white")
      par(bg = "white")


      msub <- Hmisc::subplot(hist(df$mon, breaks = c(0:12), col = ccol, xaxt = "n",
                                  bg = "white",
                                  xlab = "", main = "", cex.axis = 0.45,
                                  ylab = "Proportion",
                                  freq = FALSE, las = 1),
                             x = 0.00,  y = max(df$amax), hadj = 0, vadj = 1,
                             size = c(1.5, 0.75))
      op <- par(no.readonly=TRUE)
      par(msub)
      box()
      axis(1, at = 0.5:11.5, labels = mnths, tick = TRUE, lwd = -1,
           line = -1, cex.axis = 0.6)
      par(op)
    }
}
 #########

###################################################### end POT plot



######################################################## polar plot
if(polar){

oldpar <- par(no.readonly = TRUE)

par(mar = c(6,3,3,2))
      ptcol = "gray70"
      mcex = 0.65
      mcex1 = 0.44
      ################################################ polar plot

      par(mar = c(5,1,2,1))
      par(oma = c(1,1,1,1))

      ch_polar_plot_peaks(title = mtitle, days = df$doy, pt_col = ptcol )



      c_days <- df$doy/df$days *360
      rad_days <- deg2rad(c_days)

      cos_day <- cos(rad_days)
      sin_day <- sin(rad_days)
      mags <- df$amax/maxQ

      events <- data.frame(cos_day, sin_day, mags)


      ###############################################
      ###############################################
      plotrix::radial.plot(lengths = mags,
                           radial.pos = rad_days,
                           rp.type = "s",
                           point.symbols = qpch[1],
                           point.col = qcol [1],
                           cex =  qcex[1],
                           start = 3 * pi/2,
                           clockwise = TRUE,
                           radial.lim = c(0, 1.0),
                           add = TRUE)

      plotrix::radial.plot(lengths = mags[mindex],
                           radial.pos = rad_days[mindex],
                           rp.type = "s",
                           point.symbols = qpch[2],
                           point.col = qcol [2],
                           cex =  qcex[2],
                           start = 3 * pi/2,
                           clockwise = TRUE,
                           radial.lim = c(0, 1.0),
                           add = TRUE)

      ###############################################################
      mq1 <- mags[-mindex]

#################################  Grubbs high outliers

      ghindex <- which(gtest$tout == 1)


      if(gt0$p.value <= 0.05){

        plotrix::radial.plot(lengths = mags[mindex],
                             radial.pos = rad_days[mindex],
                             rp.type = "s",
                             point.symbols = qpch[3],
                             point.col = qcol [3],
                             cex =  qcex[3],
                             start = 3 * pi/2,
                             clockwise = TRUE,
                             radial.lim = c(0, 1.0),
                             add = TRUE)

        plotrix::radial.plot(lengths = mags[ghindex],
                             radial.pos = rad_days[ghindex],
                             rp.type = "s",
                             point.symbols = qpch[3],
                             point.col = qcol [3],
                             cex =  qcex[3],
                             start = 3 * pi/2,
                             clockwise = TRUE,
                             radial.lim = c(0, 1.0),
                             add = TRUE)


      }


      #####################################  Grubbs low outliers


      if(mg$klow != 0){

        gindex <- which(df$amax < mg$LOThresh)
        plotrix::radial.plot(lengths = mags[gindex],
                             radial.pos = rad_days[gindex],
                             rp.type = "s",
                             point.symbols = qpch[4],
                             point.col = qcol [4],
                             cex =  qcex[4],
                             start = 3 * pi/2,
                             clockwise = TRUE,
                             radial.lim = c(0, 1.0),
                             add = TRUE)

      }

      rad_days <- circular::circular(rad_days, units = "radians", modulo = "asis",
                                     type = "angles", template = "none",
                                        zero = 3*pi/2, rotation = "clock")



##########
### could make the bandwidths of 15 and 50 parameters

      res15 <- circular::density.circular(rad_days, bw=15,
                                          control.circular = list(start = 3 * pi/2, clockwise = TRUE))
      res50 <- circular::density.circular(rad_days, bw=50,
                                          control.circular = list(start = 3 * pi/2, clockwise = TRUE))

      ###########################################  plotting
      ####################### add points based on rescaled magnitude (... 1.0)
      plotrix::radial.plot(lengths = mags,
                           radial.pos = rad_days,
                           rp.type = "s",
                           point.symbols = qpch[1],
                           point.col = qcol [1],
                           cex =  qcex[1],
                           start = 3 * pi/2,
                           clockwise = TRUE,
                           radial.lim = c(0, 1.0),
                           add = TRUE)


      plotrix::radial.plot(lengths = mags[mindex],
                           radial.pos = rad_days[mindex],
                           rp.type = "s",
                           point.symbols = qpch[2],
                           point.col = qcol [2],
                           cex =  qcex[2],
                           start = 3 * pi/2,
                           clockwise = TRUE,
                           radial.lim = c(0, 1.0),
                           add = TRUE)

      lines(res50, col="red", lw = 1.5)
      lines(res15, col="green", lw = 1.5)


      f_out <- which(isout == 2)

      SR_out <- f_out

      f_out0 <-rad_days[f_out]



      lth <- rep(1,length(f_out0))

      plotrix::radial.plot(lengths = lth,
                           radial.pos = f_out0,
                           rp.type = "s",
                           point.symbols = qpch[5],
                           point.col = qcol [5],
                           cex =  qcex[3],
                           start = 3 * pi/2,
                           clockwise = TRUE,
                           radial.lim = c(0, 1.0),
                           add = TRUE)



      plotrix::radial.plot(lengths = mags[mindex],
                           radial.pos = rad_days[mindex],
                           rp.type = "s",
                           point.symbols = qpch[2],
                           point.col = qcol [2],
                           cex =  qcex[2],
                           start = 3 * pi/2,
                           clockwise = TRUE,
                           radial.lim = c(0, 1.0),
                           add = TRUE)

      c_test <- mle.vonmises(rad_days, control.circular = list(start = 3 * pi/2, clockwise = TRUE))
      mkappa <- c_test$kappa

      print(paste("kappa =", round(mkappa, 4), "κ threshold =", kthresh))


      if(mkappa <= kthresh){
################################################################## fuzzy clustering
  print("kappa < κ threshold:  doing fuzzy clustering")

        qcol1 = c("darkblue",  "blue", "red")

        ########################################### fuzzy clustering

        res.upfc <- upfc(events, centers=2)  #unsupervised possibilistic fuzzy [ppclust]
        #head(res.upfc$t)


        ################# try from 2 to 5 clusters

        c1 <- 2  #Starting number of clusters
        c2 <- 5  #Final number of clusters
        indnames <- c("PC","MPC","PE","XB","Kwon", "TSS", "CL", "FS", "PBMF","FSIL","FHV", "APD")
        indvals <- matrix(ncol=length(indnames), nrow=(c2-c1+1))
        colnames(indvals) <- indnames
        rownames(indvals) <- paste0("c=", c1:c2)
        i <- 1
        for(c in c1:c2){
          resfcm <- ppclust::fcm(x=events, centers=c, nstart=3)
          indvals[i,1] <- pc(resfcm)
          indvals[i,2] <- mpc(resfcm)
          indvals[i,3] <- pe(resfcm)
          indvals[i,4] <- xb(resfcm)
          indvals[i,5] <- kwon(resfcm)
          indvals[i,6] <- tss(resfcm)
          indvals[i,7] <- cl(resfcm)
          indvals[i,8] <- fs(resfcm)
          indvals[i,9] <- pbm(resfcm)
          indvals[i,10] <- si(resfcm)$sif
          indvals[i,11] <- fhv(resfcm)
          indvals[i,12] <- apd(resfcm)
          i <- i+1
        }




        # Display the fuzzy indices in various runs of FCM
        indvals <- round(t(indvals),3)

        # Optimal number of clusters with Fuzzy Hypervolume (FHV) index
        optk <- colnames(indvals)[which.min(indvals["FHV",])]

        k <- unname(which.min(indvals["FHV",])) + 1



        ####################  use k

        res.upfc <- upfc(events, centers = k)
        res.out <- detect.outliers(res.upfc, k = k, alpha=0.01, alpha2=0.05, tsc = "m1")


        mcases <- unname(res.out$outliers1)
        fz_out <- mcases

        forecord <- which(res.out$X[,3] ==1)


        m_max <- which.max(res.out$rowSums)

        #####################  label members of clusters
        cl_col <- c( "darkgreen","blue",  "yellowgreen","cyan", "green")
        cl_num <- unname(res.upfc$cluster)

        plotrix::radial.plot(lengths = mags,
                             radial.pos = rad_days,
                             rp.type = "s",
                             point.symbols = 19,
                             point.col = cl_col[cl_num],
                             cex =  0.75,
                             start = 3 * pi/2,
                             clockwise = TRUE,
                             radial.lim = c(0, 1.0),
                             add = TRUE)

        plotrix::radial.plot(lengths = mags[mindex],
                             radial.pos = rad_days[mindex],
                             rp.type = "s",
                             point.symbols = qpch[2],
                             point.col = qcol1 [2],
                             cex =  qcex[2],
                             start = 3 * pi/2,
                             clockwise = TRUE,
                             radial.lim = c(0, 1.0),
                             add = TRUE)

        ###############  mark odetector outliers


        plotrix::radial.plot(lengths = mags[mcases],
                             radial.pos = rad_days[mcases],
                             rp.type = "s",
                             point.symbols = 19,
                             point.col = "red",
                             cex =  qcex[1],
                             start = 3 * pi/2,
                             clockwise = TRUE,
                             radial.lim = c(0, 1.0),
                             add = TRUE)


        cltext <- c(1:k, "outlier")
        legend(0.85, 1.00 ,cltext, title = "Fuzzy Cluster",
               col = c(cl_col[1:k],"red"), pch = 19, pt.cex = qcex[1],
               cex = 0.65)




################################################################# fuzzy clustering
      }

      ltext <- c("amax", "FoR",
                 "Grubbs-high","Grubbs-low", "Timing outlier")
      legend(0.85, -0.8 ,legend = ltext, col = qcol, pch = qpch, pt.cex = qcex,
             cex = 0.65)
############################# add summary

      mtext(paste("Flood of record", maxQ, "m3/s on", maxDate, sep = " "),
       side = 1, line = 0.0, cex = mcex, outer = TRUE, adj = 0., font = 2)

if(trans == "none") {
if(length(ghindex) >= 1)
  mtext(paste("# high Grubbs outliers =", length(ghindex)),
               font = 2, col ="red", side = 1, line = -2.9, cex = mcex, outer = TRUE, adj = 0.)
if(!length(ghindex) >= 1)
    mtext(paste("# high Grubbs outliers =", length(ghindex)),
          font = 2, col ="red", side = 1, line = -2.9, cex = mcex, outer = TRUE, adj = 0.)

if(mg$klow >= 1)
   mtext(paste("# low Grubbs outliers =", mg$klow),
         font = 2, col ="red",   side = 1, line = -2.3, cex = mcex, outer = TRUE, adj = 0.)
if(mg$klow >= 1)
      mtext(paste("# low Grubbs outliers =", mg$klow),
         side = 1, line = -2.3, cex = mcex, outer = TRUE, adj = 0.)
}


if(trans == "log"){
  if(length(ghindex) >= 1)
    mtext(paste("# high Grubbs log outliers =", length(ghindex)),
          font = 2, col ="red", side = 1, line = -2.9, cex = mcex, outer = TRUE, adj = 0.)
  if(!length(ghindex) >= 1)
    mtext(paste("# high Grubbs log outliers =", length(ghindex)),
          side = 1, line = -2.9, cex = mcex, outer = TRUE, adj = 0.)

  if(mg$klow >= 1)
    mtext(paste("# low Grubbs log outliers =", mg$klow),
          font = 2, col ="red",   side = 1, line = -2.3, cex = mcex, outer = TRUE, adj = 0.)
  if(!mg$klow >= 1)
    mtext(paste("# low Grubbs log outliers =", mg$klow),
          side = 1, line = -2.3, cex = mcex, outer = TRUE, adj = 0.)
}

  if(trans == "boxcox"){
        if(length(ghindex) >= 1)
          mtext(paste("# high Grubbs BoxCox outliers =", length(ghindex)),
                font = 2, col ="red", side = 1, line = -2.9, cex = mcex, outer = TRUE, adj = 0.)
        if(!length(ghindex) >= 1)
          mtext(paste("# high Grubbs BoxCox outliers =", length(ghindex)),
               side = 1, line = -2.9, cex = mcex, outer = TRUE, adj = 0.)

        if(mg$klow >= 1)
          mtext(paste("# low Grubbs BoxCox outliers =", mg$klow),
                font = 2, col ="red",   side = 1, line = -2.3, cex = mcex, outer = TRUE, adj = 0.)
        if(mg$klow >= 1)
          mtext(paste("# low Grubbs BoxCox outliers =", mg$klow),
                side = 1, line = -2.3, cex = mcex, outer = TRUE, adj = 0.)
      }



if(length(f_out) >=1)
      mtext(paste("Timing outliers= ",length(f_out), " (", round(porcout,3),") ",
                  " n = ", length (events[,1]), sep = ""),
            font = 2, col = "red",side = 1, line  = -1.7, cex = mcex, outer = TRUE, adj = 0.)
if(!length(f_out) >=1)
      mtext(paste("Timing outliers= ",length(f_out), " (", round(porcout,3),") ",
                  " n = ", length (events[,1]), sep = ""),
            side = 1, line  = -1.7, cex = mcex, outer = TRUE, adj = 0.)

if((mkappa - kthresh) <= 0 )
      mtext(paste("fuzzy clustering: kappa =", round(mkappa, 3), "threshold =", kthresh),
            font = 2, side =1, line = -1.1, cex = mcex, outer= TRUE, adj = 0.)

if(!(mkappa - kthresh) <= 0 )
        mtext(paste("no fuzzy clustering: kappa =", round(mkappa, 3), "threshold =", kthresh),
              side =1, line = -1.1, cex = mcex, outer= TRUE, adj = 0.)


 if(POT) mtext(paste("POT events  threshold =", thresh), side = 1, line = 0.60,
               cex = 0.9, adj = 1.)
      ############################################################### end of  panel 1 outliers by Sau & Rodriguez end
par(oldpar)
    }

    df$h_grubb <- gtest$tout + 1
    df$l_grubb  <- l_grubb
    df$t_outlier <- isout
    df$fcluster <- cl_num


######################################################## return details
    mresult <- list(station = stn,
               nyears = nyears,
               maxQ = maxQ,
               mindex = mindex,
               maxDate = maxDate,
               maxdoy = maxdoy,
               transform = trans,
               WaldWolfowitz = wwo,
               MannKendall = mk,
               Pettitt = pet,
               POTmodel = fitted,
               out_hgrubb = length(ghindex),
               out_lgrubb = mg$klow,
               out_timing = length(f_out),
               fz_out = mcases,
               df = df)

       return(mresult)
  }




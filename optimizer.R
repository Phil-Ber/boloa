suppressMessages(library(xcms))
suppressMessages(library(dplyr))
### This script is used to optimize the parameters for peak detection. 
### It takes three files as input with corresponding arrays concerning peak
### spread.
###

print("Running optimization... v2")

set.seed(143241)
register(bpstart(MulticoreParam(8)))
### Manually counted peaks with the corresponding file locations
## Lower peak limit
lf <- list(c(243, 248, 253, 261, 336, 345, 470, 479, 567, 600, 672),
           c(243,249,253,261,392,461,460,478,563,600,617,629,637,653,672,710),
           c(243,248,253,262,284,392,458,468,478,493,524,549,565,575,599,622,629,638,644,653,665,672,710)
)
## Upper peak limit
uf <- list(c(248, 253, 258, 268, 340, 347, 474, 482, 571, 603, 676),
           c(247,253,257,269,394,463,473,482,569,603,620,632,641,655,676,713),
           c(248,253,258,268,289,394,460,473,482,496,526,553,569,579,604,624,631,640,646,655,667,675,712)
)
## Array of paths to analyzed files.
msfile <- c("/exports/nas/berends.p/boloa/massascans/062475b654c71a3df8103a49dcf4102fc7f503d1c886f8d6a91371ce.mzXML",
            '/exports/nas/berends.p/boloa/massascans/0f4a154d2423caf8967d1d0275438c6c0d59c51e3d5fb5cec90b6cf4.mzXML',
            "/exports/nas/berends.p/boloa/massascans/107c76e6b73bbb1c899aaac73696af6c24fbafa2670c152827a81cc8.mzXML"
)
## How wide a peak can be
peak_margin <<- 4

pd <- data.frame(sample_name = c("Blanco2_1", "T0128_1150", 'QC4_6'), sample_group = c("Blanco2", "A", "QC4"), stringsAsFactors = FALSE)
raw_data <- readMSData(files = msfile, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")
raw_data <<- filterRt(raw_data, c(min(rtime(raw_data)), 750))

detection <- function(detection_algo, params) {
  tryCatch(
    {
      if (detection_algo == 1) {
        xset <<- suppressMessages(findChromPeaks(raw_data, param = CentWaveParam(
          ppm = params[1],
          peakwidth = c(params[2], params[3]),
          snthresh = params[4],
          prefilter = c(round(params[5]), params[6]),
          mzCenterFun = c("wMean", "mean", "apex", "wMeanApex3", "meanApex3")[round(params[7])],
          integrate = 1L,
          mzdiff = params[8],
          fitgauss = as.logical(round(params[9])),
          noise = params[10],
          verboseColumns = as.logical(round(params[11])
          ))))
      } else if (detection_algo  == 2) {
        xset <<- findChromPeaks(raw_data, param = MatchedFilterParam(params))
      } else {
        xset <<- findChromPeaks(raw_data, param = MassifquantParam(params))
      }
      if (is.null(nrow(chromPeaks(xset)))) {
        return(as.data.frame(t(c(1, 1, 1, 1, params))))
      }
      fres <- c()
      for (file_index in 1:length(msfile)) {
        n_corr <- 0
        peaksf <- chromPeaks(xset, bySample = 1)[file_index]
        for (ci in 1:length(lf[[file_index]])) {
          # Check if all the peaks are solely within the given regions.
          corr_frame <- peaksf[[1]][between(peaksf[[1]][,"rtmin"], lf[[file_index]][ci] - peak_margin, uf[[file_index]][ci] + peak_margin) & between(peaksf[[1]][,"rtmax"], lf[[file_index]][ci] - peak_margin, uf[[file_index]][ci] + peak_margin),]
          cor_det <- replace(length(corr_frame) / length(corr_frame), is.na(length(corr_frame) / length(corr_frame)), 1)
          n_corr <- n_corr + cor_det # cor_det is a binary value, either 0 or 1, 1 signifies a correct region is detected.
        }
        n_corr <- n_corr - (-nrow(peaksf[[1]]) + n_corr)**2
        fres <- c(fres, abs(1-length(lf[[file_index]])/(abs(n_corr)+length(lf[[file_index]]))-0.5))
      }
      return(as.data.frame(t(c(mean(fres), fres, params))))
    },
    error = function(cnd){
      return(as.data.frame(t(c(1, 1, 1, 1, params))))
    }
  )
}


start_vals <- c(93, 5, 25, 100, 1, 300, 1, 100, 1, 1000, 0)
options(warn=-1)
edf <- detection(1, start_vals)
for (i in 1:1000) {
  print(i)
  minpeak <- runif(1,1.1,100)
  randvals <- c(
    runif(1, 1.1, 200), #ppm
    minpeak, #min_peakwidth
    runif(1, minpeak, 101+minpeak), #min+val = max_peakwith
    runif(1, 1.1, 100), #snthresh
    sample(c(1,2,3), 1), #prefilter
    runif(1, 1.1, 400), #val_of_prefilter
    sample(c(1,2,3,4,5), 1), #mzCenterFun
    runif(1, -1.1, 1), #mzdiff
    sample(c(0,1), 1), #fitgauss
    runif(1, 0.1, 10000), #noise
    sample(c(0,1), 1) #verbosecolumns
  )
  edf <- rbind(edf, detection(1, randvals))
  #edf <- edf[edf[,1] != 1,]
}
write.csv(edf, "/exports/nas/berends.p/boloa/optimize_screen_4.csv", row.names=FALSE)

tries <- 0
prev_add <- data.frame(c(2,2,3))
while (mean(edf[,1]) >= 0.10 | tries < 100) {
  edf <- tail(edf, n = 1000) #For testing
  end_vals <- c(t(colMeans(edf[,5:length(edf)])))
  end_vals <- (end_vals + c(t(edf[order(edf[,1]),][1,5:15]))) / 2
  #message("N TRIES:")
  #message(tries)
  if (tries %% 5 == 0) {
    npars <- c()
    for (i in 1:length(end_vals)) {
      GSfactor <- -cor(edf[,i+4], edf[,1])
      npars <- c(npars, end_vals[i] * (GSfactor**2 + 1))
    }
    
  }
  else if (tries %% 12 == 0) {
    npars <- c()
    for (i in 1:length(end_vals)) {
      GSfactor <- -cor(edf[,i+4], edf[,1])
      npars <- c(npars, end_vals[i] * (GSfactor**3 + 1))
    }
    edf <- edf[edf[,1] != 1,]
  }
  else if (tries %% 49 == 0) {
    npars <- c()
    for (i in 1:length(end_vals)) {
      GSfactor <- -cor(edf[,i+4], edf[,1])
      npars <- c(npars, end_vals[i] * (GSfactor**4 + 1))
    }
    edf <- edf[edf[,1] != 1,]
  }
  else {
    npars <- c()
    for (i in 1:length(end_vals)) {
      GSfactor <- -cor(edf[,i+4], edf[,1])
      npars <- c(npars, end_vals[i] * (GSfactor + 1))
    }
  }
  new_add <- detection(1, npars)
  #message(new_add)
  if (prev_add[,1][1] == new_add[,1]) {
    tries <- tries + 1
  } else {
    tries <- 0 
  }
  edf <- rbind(edf, new_add)
  prev_add <- new_add
  write.table(new_add, file = "/exports/nas/berends.p/boloa/optimize_screen_4.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
}
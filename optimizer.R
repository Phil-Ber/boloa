#!/usr/bin/env Rscript
suppressMessages(library(xcms))
suppressMessages(library(dplyr))
### This script is used to optimize the parameters for peak detection. 
### It takes three files as input with corresponding arrays concerning peak
### spread.
###

print("Running optimization... v2")

set.seed(2023)
register(bpstart(MulticoreParam(8)))
### Manually counted peaks with the corresponding file locations
# Lower peak limit
lf <- list(c(261, 336, 345, 470, 479, 567, 600, 672), #F1
           c(261,392,461,460,478,563,600,617,629,637,653,672,710), #F2
           c(262,284,392,458,468,478,493,524,549,565,575,599,622,629,638,644,653,665,672,710) #F3
)
## Upper peak limit
uf <- list(c(268, 340, 347, 474, 482, 571, 603, 676), #F1
           c(269,394,463,473,482,569,603,620,632,641,655,676,713), #F2
           c(268,289,394,460,473,482,496,526,553,569,579,604,624,631,640,646,655,667,675,712) #F13
)
## Array of paths to analyzed files.
msfile <- c("/massascans/062475b654c71a3df8103a49dcf4102fc7f503d1c886f8d6a91371ce.mzXML", #F1
            '/massascans/0f4a154d2423caf8967d1d0275438c6c0d59c51e3d5fb5cec90b6cf4.mzXML', #F2
            "/massascans/107c76e6b73bbb1c899aaac73696af6c24fbafa2670c152827a81cc8.mzXML" #F3
)

peak_margin <<- 1

pd <- data.frame(sample_name = c("Blanco2_1", "T0128_1150", 'QC4_6'), sample_group = c("Blanco2", "A", "QC4"), stringsAsFactors = FALSE)
raw_data <- readMSData(files = msfile, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")
raw_data <<- filterRt(raw_data, c(260, 750))

detection <- function(detection_algo, params) {
  # This function performs peak detection. It requires two variables:
  #   detection algo, which is the algorithm which has been used for peak detection
  #   prarams, the parameters for the algorithm.
  #
  # A dataframe is returned with three scores and the used parameters
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
          cor_det <- replace(length(corr_frame) / length(corr_frame), is.na(length(corr_frame) / length(corr_frame)), 0)
          n_corr <- n_corr + cor_det # cor_det is a binary value, either 0 or 1, 1 signifies a correct region is detected.
        }
        n_corr <- n_corr - (nrow(peaksf[[1]]) + n_corr)
        fres <- c(fres, pmax(0, 1 - abs(n_corr + length(lf[[file_index]])) / (n_corr - length(lf[[file_index]]))^2))
        fres[is.nan(fres)] <- 1
      }
      return(as.data.frame(t(c(mean(fres), fres, params))))
    },
    error = function(cnd){
      return(as.data.frame(t(c(1, 1, 1, 1, params))))
    }
  )
}

# Start values (not really necessary)
start_vals <- c(93, 5, 25, 100, 1, 300, 1, 100, 1, 1000, 0)
# start_vals <- c(93, 25, 100, 1, 300, 1, 100, 1, 1000, 0)
options(warn=-1)
edf <- detection(1, start_vals)
# Generate 1000 parameter and score combinations
for (i in 1:1000) {
  print(i)
  print(Sys.time())
  minpeak <- runif(1,1.1,10)
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
  edf <- edf[edf[,1] != 1,]
}

# Save to file for further usage
write.csv(edf, "randomgen_param_for_optimization.csv", row.names=FALSE)

# Load randomly generated parameters
edf <- read.csv("andomgen_param_for_optimization.csv", header
                =TRUE, stringsAsFactors=FALSE)

# Remove NAs
rng <- na.omit(edf)
edf <- na.omit(edf)
# Unused best parameters
best_100 <- head(edf[order(edf[,1]),], n = 100)
tries <- 0
edf <- edf[edf[,1] < 1,]
prev_add <- head(best_100, n = 1)
current_best <- c(t(edf[order(edf[,1]),][1,]))

# Optimize parameters loop
while (mean(edf[,1]) >= 0.10 | tries < 100) {
  edf <- na.omit(edf)
  previous_params <- c(t(prev_add))
  
  # Replace best if previous parameters were better than the previous best
  if (previous_params[1] < current_best[1]) {
    current_best <- previous_params
  }
  
  # Take 1000 last parameters
  edf <- tail(edf, n = 1000)
  mean_over_edf <- c(t(colMeans(edf[,5:length(edf)])))
  #end_vals <- (mean_over_edf + current_best[5:15] + previous_params[5:15]) / 3
  # Calculate new parameters from mean between average parameters and best
  end_vals <- (mean_over_edf + c(t(edf[order(edf[,1]),][1,5:15]))) / 2
  
  # Calculate correlation factors
  GSfactor <- -(cor(edf, edf[,1])[5:ncol(edf)])
  
  # Apply factors to parameters
  if (tries %% 2 == 0) {
    npars <- end_vals * (GSfactor**2 + 1)
  }
  else if (tries %% 5 == 0) {
    npars <- end_vals * (GSfactor**3 + 1)
  }
  else {
    npars <- end_vals * (GSfactor + 1)
  }
  new_add <- detection(1, npars)
  #message(new_add)
  if (prev_add[,1][1] == new_add[,1]) {
    tries <- tries + 1
  } else {
    #tries <- 0 
  }
  # Add to dataframe of scores and parameters
  edf <- rbind(edf, new_add)
  prev_add <- new_add
  write.table(new_add, file = "optimize_screen_17.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
}
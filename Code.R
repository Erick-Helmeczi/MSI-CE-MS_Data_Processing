library("xcms")
library("doParallel")
library("BiocParallel")
library("tidyverse")
library("pracma")
library("dplyr")
library("DescTools")

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

file <- "20210425_40mM.mzML"

registerDoParallel(3) 
register(DoparParam(), default = TRUE)

rawdata <- readMSData(
  file = file,
  pdata = NULL,
  msLevel = 1,
  verbose = isMSnbaseVerbose(),
  centroided. = FALSE,
  smoothed. = FALSE,
  cache. = 0,
  mode =  "onDisk"
)

min = 509.3738 - (509.3738/20000)
max = 509.3738 + (509.3738/20000)

mzrange <- cbind(min, max)
rtrange <- cbind(120, 600)

EIE <- chromatogram(rawdata, 
                    mz = mzrange, 
                    rt = rtrange,
                    aggregationFun = "sum",
                    missing = NA_real_,
                    msLevel = 1L)


rtime_Vector <- unname(EIE[1]@rtime/60)
Intensity_Vector <- unname(EIE[1]@intensity)

data_frame <- data.frame(rtime_Vector, Intensity_Vector)
tdata_frame <- subset(data_frame, rtime_Vector >= min)
tdata_frame <- subset(tdata_frame, rtime_Vector <= max)

Peaks <- findpeaks(Intensity_Vector,
                   nups  = 3,
                   ndowns = 3,
                   zero = 0,
                   minpeakheight = 15000,
                   minpeakdistance = 6,
                   npeaks = 6,
                   sortstr = FALSE
)

Peaks

Min_MT_Vec <- data_frame[Peaks[,3],1]
Max_MT_Vec <- data_frame[Peaks[,4],1]
Apex_MT_Vec <- data_frame[Peaks[,2],1]

Peak_df <- data.frame(Min = Min_MT_Vec, Max = Max_MT_Vec, Apex = Apex_MT_Vec)

Peak_df

Total_Area <- AUC(rtime_Vector, 
                  Intensity_Vector, 
                  from = Peak_df[1,1], 
                  to = Peak_df[1,2], 
                  method = "trapezoid")

Baseline <- (Peak_df[1,2] - Peak_df[1,1])*data_frame[Peaks[1,3],2]

Area <- Total_Area - Baseline

ggplot(data = data_frame, aes(x = rtime, y = intensity)) +
  geom_line() +
  labs(title = "EIE m/z = 509.3738",
       x = "Migration Time (Minutes)",
       y = "Ion Counts") +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        strip.text = element_text(face = "italic"),
        text = element_text(size = 16)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_area(data = tdata_frame, alpha = 0.5, fill = "blue") +
  theme(legend.position = "none")
library("xcms")
library("doParallel")
library("BiocParallel")
library("tidyverse")
library("pracma")
library("dplyr")
library("DescTools")
library("geiger")
library("MESS")
library("mhsmm")
library("xlsx")
library("openxlsx")

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

file <- "20210425_10mM_60C.mzML"

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

Mass_df <- read.xlsx("Mass List.xlsx")

for (i in 1:nrow(Mass_df)-1){

min = Mass_df[i+1,2] - (Mass_df[i+1,2]/20000)
max = Mass_df[i+1,2] + (Mass_df[i+1,2]/20000)

mzrange <- cbind(min, max)
rtrange <- cbind(100, 840)

EIE <- chromatogram(rawdata, 
                    mz = mzrange, 
                    rt = rtrange,
                    aggregationFun = "sum",
                    missing = NA_real_,
                    msLevel = 1L)


rtime_Vector <- unname(EIE[1]@rtime)
Intensity_Vector <- unname(EIE[1]@intensity)

data_frame <- data.frame(rtime_Vector, Intensity_Vector)

Smooth <- with(data_frame, 
            ksmooth(rtime_Vector, Intensity_Vector, kernel = "box", bandwidth = 10))

data_frame <- data.frame(rtime_Vector = Smooth$x, Intensity_Vector = Smooth$y)
rtime_Vector <- data_frame$rtime_Vector
Intensity_Vector <- data_frame$Intensity_Vector

Peaks <- findpeaks(Intensity_Vector,
                   nups  = 3,
                   ndowns = 3,
                   zero = "+",
                   minpeakheight = 15000,
                   minpeakdistance = 6,
                   npeaks = 6,
                   sortstr = FALSE
)

Min_MT_Vec <- data_frame[Peaks[,3],1]
Max_MT_Vec <- data_frame[Peaks[,4],1]
Apex_MT_Vec <- data_frame[Peaks[,2],1]

Peak_df <- data.frame(Min = Min_MT_Vec, Max = Max_MT_Vec, Apex = Apex_MT_Vec)
Peak_df <- Peak_df[order(Peak_df$Min),]
Peak_df

Results_df <- data.frame(Area = "")

for (i in 1:nrow(Peak_df)){
  Area <- auc(rtime_Vector, 
              Intensity_Vector, 
              from = Peak_df[i,1],
              to = Peak_df[i,2],
              type = "linear",
              absolutearea = TRUE)
  Results_df[i,1] <- Area
}

Results_df <- cbind(Results_df, Migration_Time = Peak_df[,3]/60)

Name <- Mass_df[i+1,1]

geom_area_df <- data.frame()
for (i in 1:nrow(Peak_df)){
  geom_area_df_temp <- data_frame %>% filter(rtime_Vector >= Peak_df[i,1], rtime_Vector <= Peak_df[i,2])
  geom_area_df <- rbind(geom_area_df, geom_area_df_temp)
}

ggplot(data = data_frame, aes(x = rtime_Vector/60, y = Intensity_Vector)) +
  geom_line() +
  labs(title = Name,
       x = "Migration Time (Minutes)",
       y = "Ion Counts") +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        strip.text = element_text(face = "italic"),
        text = element_text(size = 16)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_area(data = geom_area_df, alpha = 0.5, fill = "blue") +
  theme(legend.position = "none")

ggsave(filename = Name.pdf,
       plot = last_plot())

write.xlsx(Results_df,
           file = "MSI-CE-MS_Data_Processing.xlsx",
           sheetName = Name,
           append = TRUE)

}

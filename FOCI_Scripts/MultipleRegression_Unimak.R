#read data into variable
pathname <- "/Users/bell/in_and_outbox/2016/stabeno/feb/unimakwinds_narr/"
filename <- "NARR_Winds_shumigans_UnimakP_tran1996.csv"
datavar <- read.csv(paste(pathname,filename,sep = ""))

#attach data variable
attach(datavar)

#two predictor model
twoPredictorModel <- lm(u_curr_cms ~ WU_422 + WV_423, datavar)

# summary
twoPredictorModel
summary(twoPredictorModel)

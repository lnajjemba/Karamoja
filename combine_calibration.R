library(dplyr)

outputs <- NULL #initiate blank output object
js.set <- 1:100 #number of files (100, or 5000). If you just ran calibration on your computer (not the server), this is 1.

#load output from each file
for(js in js.set){ 
	tres <- read.csv(paste0('calibration/summary_output_', js, '.csv'))
	outputs <- rbind(outputs,tres)
}

#read in calibration targets
targets <- read.csv('calibration_targets.csv', row.names=1)

#Option 1: take results if they are within 25% of the means of all targets
threshold <- 0.25 
select1 <- which(outputs$prev > targets["mean", "prev"]*(1-threshold) &
                   outputs$prev < targets["mean", "prev"]*(1+threshold) & 
                   outputs$mort > targets["mean", "mort"]*(1-threshold) &
                   outputs$mort < targets["mean", "mort"]*(1+threshold) &
                   outputs$notif_2017 > targets["mean", "notif_2017"]*(1-threshold) &
                   outputs$notif_2017 < targets["mean", "notif_2017"]*(1+threshold) &
                   outputs$notif_2022 > targets["mean", "notif_2022"]*(1-threshold) &
                   outputs$notif_2022 < targets["mean", "notif_2022"]*(1+threshold) &
                   outputs$prop_notif_u15 > targets["mean", "prop_notif_u15"]*(1-threshold) &
                   outputs$prop_notif_u15 < targets["mean", "prop_notif_u15"]*(1+threshold) &
                   outputs$rel_inc_undernut > targets["mean", "rel_inc_undernut"]*(1-threshold) &
                   outputs$rel_inc_undernut < targets["mean", "rel_inc_undernut"]*(1+threshold) &
                   outputs$prop_prev_sub > targets["mean", "prop_prev_sub"]*(1-threshold) &
                   outputs$prop_prev_sub < targets["mean", "prop_prev_sub"]*(1+threshold) &
                   outputs$inc > targets["mean", "inc"]*(1-threshold) &
                   outputs$inc < targets["mean", "inc"]*(1+threshold)
)



#Option 2: take results if they are within CIs of all targets
select2 <- which(outputs$prev > targets["low", "prev"] &
                   outputs$prev < targets["high", "prev"] & 
                   outputs$mort > targets["low", "mort"] &
                   outputs$mort < targets["high", "mort"] &
                   outputs$notif_2017 > targets["low", "notif_2017"] &
                   outputs$notif_2017 < targets["high", "notif_2017"] &
                   outputs$notif_2022 > targets["low", "notif_2022"] &
                   outputs$notif_2022 < targets["high", "notif_2022"] &
                   outputs$prop_notif_u15 > targets["low", "prop_notif_u15"] &
                   outputs$prop_notif_u15 < targets["high", "prop_notif_u15"] &
                   outputs$rel_inc_undernut > targets["low", "rel_inc_undernut"] &
                   outputs$rel_inc_undernut < targets["high", "rel_inc_undernut"] &
                   outputs$prop_prev_sub > targets["low", "prop_prev_sub"] &
                   outputs$prop_prev_sub < targets["high", "prop_prev_sub"] &
                   outputs$inc > targets["low", "inc"] &
                   outputs$inc < targets["high", "inc"]
)

#Option 3: weight results based on likelihood and resample, with replacement
#this will require revisions once targets are finalized
log_like_prev <- dnorm(x=outputs$prev, mean=targets["mean", "prev"], 
                   sd=(targets["high", "prev"] - targets["low", "prev"])/(2*1.96), log=T)
log_like_mort <- dnorm(x=outputs$mort, mean=targets["mean", "mort"], 
                   sd=(targets["high", "mort"] - targets["low", "mort"])/(2*1.96), log=T)
log_like_notif_2017 <- dnorm(x=outputs$notif_2017, mean=targets["mean", "notif_2017"], 
                             sd=(targets["high", "notif_2017"] - targets["low", "notif_2017"])/(2*1.96), log=T)
log_like_notif_2022 <- dnorm(x=outputs$notif_2022, mean=targets["mean", "notif_2022"], 
                   sd=(targets["high", "notif_2022"] - targets["low", "notif_2022"])/(2*1.96), log=T)
log_like_prop_notif_u15 <- dnorm(x=outputs$prop_notif_u15, mean=targets["mean", "prop_notif_u15"], 
                             sd=(targets["high", "prop_notif_u15"] - targets["low", "prop_notif_u15"])/(2*1.96), log=T)
log_like_rel_inc_undernut <- dnorm(x=outputs$rel_inc_undernut, mean=targets["mean", "rel_inc_undernut"], 
                    sd=(targets["high", "rel_inc_undernut"] - targets["low", "rel_inc_undernut"])/(2*1.96), log=T)
log_like_prop_prev_sub <- dnorm(x=outputs$prop_prev_sub, mean=targets["mean", "prop_prev_sub"], 
                    sd=(targets["high", "prop_prev_sub"] - targets["low", "prop_prev_sub"])/(2*1.96), log=T)
log_like_inc <- dnorm(x=outputs$inc, mean=targets["mean", "inc"], 
                       sd=(targets["high", "inc"] - targets["low", "inc"])/(2*1.96), log=T)
log_like_all <- log_like_prev + log_like_mort + log_like_notif_2022 + log_like_notif_2017 +
  log_like_prop_notif_u15 + log_like_rel_inc_undernut + log_like_prop_prev_sub + log_like_inc
like_all <- exp(log_like_all)
weights <- like_all/sum(like_all)
outputs <- cbind("id"=1:length(weights), outputs, "likelihood"=like_all, "weight"=weights)
select3 <- sample(outputs$id, size=1000, prob=weights, replace=T)

#collapse to version with only unique samples, but keep track of how many times we'll need to duplicate each sample
outputs3 <- outputs[select3,]
outputs3_sub <- outputs3 %>% group_by_all() %>% summarise(number=n())
select3_sub <- unique(outputs3$id)

#save to file - also add weights of 1 for each parameter set to outputs1 and outputs2
write.csv(select1, file='calibration/select_threshold_25.csv') #indices 
write.csv(select2, file='calibration/select_in_bounds.csv') #indices 
write.csv(select3, file='calibration/select_weights.csv') #indices 
write.csv(outputs[select1,] %>% mutate(number=1), file='calibration/select_results_threshold_25.csv', row.names=F) #actual parameters and outputs
write.csv(outputs[select2,] %>% mutate(number=1), file='calibration/select_results_in_bounds.csv', row.names=F) #actual parameters and outputs
write.csv(outputs3_sub, file='calibration/select_results_weights.csv', row.names=F) #actual parameters and outputs


select1.pop <- vector('list', length(select1))
sims_per_file <- nrow(tres) #total number of parameter sets simulated in each calibration run
file_old <- -1
if(length(select1)!=0) {
  for(js in 1:length(select1)){
    file <- ceiling(select1[js]/sims_per_file) 
    index <- select1[js] - sims_per_file*(file-1)
    if(file_old!=file) {
      load(paste0('calibration/model_output_' , file, '.rda'))
    }
    select1.pop[[js]] <- cbind("id"=select1[js], result[[index]]$y)
    file_old <- file
  }
}
save(select1.pop, file='calibration/select_pop_threshold_25.rda')

select2.pop <- vector('list', length(select2))
sims_per_file <- nrow(tres) #total number of parameter sets simulated in each calibration run
file_old <- -1
if(length(select2)!=0) {
  for(js in 1:length(select2)){
    file <- ceiling(select2[js]/sims_per_file) 
    index <- select2[js] - sims_per_file*(file-1)
    if(file_old!=file) {
      load(paste0('calibration/model_output_' , file, '.rda'))
    }
    select2.pop[[js]] <- cbind("id"=select2[js], result[[index]]$y)
    file_old <- file
  }
}

save(select2.pop, file='calibration/select_pop_in_bounds.rda')

select3.pop <- vector('list', length(select3_sub))
sims_per_file <- nrow(tres) #total number of parameter sets simulated in each calibration run
file_old <- -1
for(js in 1:length(select3_sub)){
  file <- ceiling(select3_sub[js]/sims_per_file) 
  index <- select3_sub[js] - sims_per_file*(file-1)
  if(file_old!=file) {
    load(paste0('calibration/model_output_' , file, '.rda'))
  }
  select3.pop[[js]] <- cbind("id"=select3_sub[js], result[[index]]$y)
  file_old <- file
}
save(select3.pop, file='calibration/select_pop_weights.rda')




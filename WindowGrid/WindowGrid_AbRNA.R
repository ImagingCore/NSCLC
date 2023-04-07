
# LOAD LIBRARIES
library(tidyverse)
library(data.table)

##### LOAD DATASET AND PREP #####

#registered
filepath = c('/Volumes/LaCie/Data/IO_Lab/Registered Images/A1_ObjectTables/')
savehere <- '/Volumes/LaCie/Data/IO_Lab/Registered Images/A2_GridWindows/'

if (!dir.exists(savehere)){
  dir.create(path=savehere)
}

# ----- NEIGHBORHOOD DISTANCE -----
window_size = 50 # um

# ----- Load Files and Wrangle -----

# Grab all filenames of files to process
files <- list.files(path = filepath, pattern = '.csv')

run_start <- Sys.time()

for (file in files) {
  
  Slide.id.Ab <- str_extract(file, pattern = '\\de\\d\\d_\\d\\d') 
  
  dt <- fread(paste0(filepath,file))
  dt <- setDT(dt)
  
  setkey(dt, 'X','Y')
  
  # registered
  markers <- c('CXCL13', 'CXCL10.11', 'CD3E', 'IFNG','CK.ish', 'NegCells', 'CD8', 'Ki67', 'PDL1', 'PD1', 'TCF7', 'CK.ab')
  
  window_displacement = window_size # 50 microns
  
  # number of windows
  n_x <- round((max(dt$X) - min(dt$X)) / (window_displacement) )
  n_y <- round((max(dt$Y) - min(dt$Y)) / (window_displacement) )
  n_win = n_x * n_y
  print(paste0('Number of windows:  ', n_win))
  
  start_time <- Sys.time()
  
  map_cellXY_to_window_center50 <- function(current_coordinate, first_target){
    lastTwoDigits <- current_coordinate %% 100
    first_border <- round((first_target-25) %% 100, 3)
    second_border <- round((first_target+25) %% 100, 3)
    second_target <- round((first_target + 50) %% 100, 3)
    
    if( first_border >= 50){
      if( first_target > first_border){
        if(first_border <= lastTwoDigits){
          center_assignment <- current_coordinate - lastTwoDigits + first_target
        } else if (lastTwoDigits <= second_border){
          center_assignment <- current_coordinate - lastTwoDigits - 100 + first_target
        } else { # second_border < current_coordinate < first_border
          center_assignment <- current_coordinate - lastTwoDigits + second_target
        }
      } else { # first target < first border
        if(first_border <= lastTwoDigits){
          center_assignment <- current_coordinate - lastTwoDigits + 100 +  first_target
        } else if (lastTwoDigits <= second_border){
          center_assignment <- current_coordinate - lastTwoDigits + first_target
        } else { # second_border < current_coordinate < first_border 
          center_assignment <- current_coordinate - lastTwoDigits + second_target
        }
      }  
    } else{ # first_border < 50
      if(second_target >= second_border){
        if(first_border >= lastTwoDigits){
          center_assignment <- current_coordinate - lastTwoDigits - 100 + second_target
        } else if (lastTwoDigits >= second_border){
          center_assignment <- current_coordinate - lastTwoDigits + second_target
        }
        else { # first_border < current_coordinate < second_border
          center_assignment <- current_coordinate - lastTwoDigits + first_target
        }
      } else { # second target < second border
        if(first_border >= lastTwoDigits){
          center_assignment <- current_coordinate - lastTwoDigits + second_target
        } else if (lastTwoDigits >= second_border){
          center_assignment <- current_coordinate - lastTwoDigits + 100 + second_target
        }
        else { # first_border < current_coordinate < second_border
          center_assignment <- current_coordinate - lastTwoDigits + first_target
        }
      }
    } # // Else (done)
    return(center_assignment)
  } # // map 50
  
  allcell_maximum_X <- max(dt$X)
  allcell_maximum_Y <- max(dt$Y) 
  target_x <- round((min(dt$X)+50) %% 100, 3)
  target_y <- round((min(dt$Y)+50) %% 100, 3)
  
  # apply function over all X and Y coordinates
  all_centerX <- lapply(dt$X, map_cellXY_to_window_center50, first_target = target_x)
  all_centerY <- lapply(dt$Y, map_cellXY_to_window_center50, first_target = target_y)
  
  dt[ , winX := mapply(all_centerX, FUN = function(x) x)][, winY := mapply(all_centerY, FUN = function(y) y)]
  
  print("All cells labelled. Binding four tables ... ")
  
  print("Summarizing ... ")
  
  nhood_counts <- dt %>% group_by(winX, winY) %>% mutate('Counts All Nhood Cells' = n()) %>%  select('Counts All Nhood Cells') %>% summarise_all(max)
  
  dt <- dt %>% 
    
    # registered
    select(CD8, CK.ab, TCF7, Ki67, PDL1, PD1, CK.ish, CXCL13, IFNG, CXCL10.11, CD3E, winX, winY) 
  
  dt[CD8==1 & PD1==1,  `Counts CD8+PD1+` := 1]
  dt[CD8==1 & Ki67==1, `Counts CD8+Ki67+`:= 1]
  dt[CD8==1 & TCF7==1, `Counts CD8+TCF7+`:= 1]
  
  dt[CD8==1 & TCF7==1 & Ki67==1, `Counts CD8+TCF7+Ki67+` := 1]
  dt[CD8==1 & PD1==1 & Ki67==1,  `Counts CD8+PD1+Ki67+`  := 1]
  dt[CD8==1 & PD1==1 & TCF7==1,  `Counts CD8+PD1+TCF7+`  := 1]
  
  dt[CD8==1 &PD1==1 & Ki67==1 & TCF7==1, `Counts CD8+PD1+Ki67+TCF7+` := 1]
  dt[CD8==1 &PD1==1 & Ki67==1 & TCF7==0, `Counts CD8+PD1+Ki67+TCF7-` := 1]
  dt[CD8==1 &PD1==1 & Ki67==0 & TCF7==1, `Counts CD8+PD1+Ki67-TCF7+` := 1]
  dt[CD8==1 &PD1==1 & Ki67==0 & TCF7==0, `Counts CD8+PD1+Ki67-TCF7-` := 1]
  dt[CD8==1 &PD1==0 & Ki67==1 & TCF7==1, `Counts CD8+PD1-Ki67+TCF7+` := 1]
  dt[CD8==1 &PD1==0 & Ki67==1 & TCF7==0, `Counts CD8+PD1-Ki67+TCF7-` := 1]
  dt[CD8==1 &PD1==0 & Ki67==0 & TCF7==1, `Counts CD8+PD1-Ki67-TCF7+` := 1]
  dt[CD8==1 &PD1==0 & Ki67==0 & TCF7==0, `Counts CD8+PD1-Ki67-TCF7-` := 1]
  
  dt[CD8==0 & PD1==1, `Counts CD8-PD1+`  := 1]
  dt[CD8==0 & TCF7==1, `Counts CD8-TCF7+` := 1]
  dt[CK.ab==0 & TCF7==1, `Counts CK.ab-TCF7+` := 1]
  
  dt[CK.ab == 1 & PDL1==1, `Counts CK.ab+PDL1+`:= 1]
  dt[CK.ab == 1 & PDL1==0, `Counts CK.ab+PDL1-`:= 1] 
  dt[CK.ab == 0 & PDL1==1, `Counts CK.ab-PDL1+`:= 1]
  dt[CK.ab == 0 & PDL1==0, `Counts CK.ab-PDL1-`:= 1] 
  dt[CK.ab == 1 & TCF7==1, `Counts CK.ab+TCF7+`:= 1]
  
  # RNA
  dt[CK.ish == 1 & CXCL10.11==1, `Counts CK.ish+CXCL10.11+`:= 1]
  dt[CK.ish == 1 & CXCL10.11==0, `Counts CK.ish+CXCL10.11-`:= 1]
  dt[CK.ish == 0 & CXCL10.11==1, `Counts CK.ish-CXCL10.11+`:= 1]
  dt[CK.ish == 0 & CXCL10.11==0, `Counts CK.ish-CXCL10.11-`:= 1]
  dt[CK.ish == 1 & CXCL13==1, `Counts CK.ish+CXCL13+`   := 1]
  
  dt[CD3E==1 & CXCL13==1, `Counts CD3E+CXCL13+`     := 1]
  dt[CD3E==1 & CXCL10.11==1, `Counts CD3E+CXCL10.11+`  := 1]
  dt[CD3E==1 & IFNG==1 & CXCL13==1, `Counts CD3E+IFNG+CXCL13+`:= 1]
  
  # registered
  setnames(dt, c('CD8', 'CK.ab', 'TCF7', 'Ki67', 'PDL1', 'PD1', 'CK.ish', 'CXCL13', 'IFNG', 'CXCL10.11', 'CD3E'),
           c('Counts CD8+', 'Counts CK.ab+', 'Counts TCF7+', 'Counts Ki67+', 'Counts PDL1+', 'Counts PD1+',
             'Counts CK.ish+', 'Counts CXCL13+', 'Counts IFNG+', 'Counts CXCL10.11+', 'Counts CD3E+') )
  
  dt[is.na(dt)] <- 0
  dt <- dt %>% group_by(winX, winY) %>% summarise_all(sum) %>% ungroup()
  
  completed_table <- full_join(dt, nhood_counts) 
  
  completed_table <- as.data.table(completed_table)
  completed_table[is.na(completed_table)] <- 0
  
  #----------------------
  end_time <- Sys.time()
  print(end_time - start_time)
  
  completed_table <- completed_table %>% add_column(`Slide.id.Ab` = Slide.id.Ab) %>% 
    mutate(Window.ID = row_number()) %>% filter(winX < allcell_maximum_X & winY < allcell_maximum_Y)
  
  setwd(savehere)
  fwrite(completed_table, file = paste0('dt_', Slide.id.Ab,  '.csv'), sep=',')   
  
  rm(completed_table)
  rm(dt)
  rm(nhood_counts)
  
} # I (FILE) // patient loop

run_end <- Sys.time()
print(paste0("The whole run of ", length(files), " took ", end_time - start_time, " minutes."))



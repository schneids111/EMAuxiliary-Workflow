library(lme4)
library(rblimp)
library(readr)
source("https://raw.githubusercontent.com/schneids111/EMAuxiliary-Workflow/main/EMAuxiliary/EMAuxiliary.R")

#### EDIT THIS LINE: Set path to your local Blimp installation
blimp_path <- "C:/Program Files/Blimp/blimp.exe"
rblimp::set_blimp(blimp_path)
#### EDIT THIS LINE: Folder where your data are stored
data_path <- "C:/path/to/your/data"
#### EDIT THIS LINE: Folder where outputs should be saved
output_path <- "C:/path/to/save/outputs"

data1 <- read_csv(file.path(data_path, "study_data.csv"))

############################################################
###multilevel model without any auxiliary variables
############################################################
EMAux_noaux <- EMAuxiliary(
  negaff ~  physact1*physact1.mean + (1 + physact1 | idn),
  data = data1,
  id = 'idn',
  center_group = 'physact1',
  center_grand = 'physact1.mean', 
  seed = 90291,
  burn = 6000,
  iter = 10000,
  chains = 4
)
cat(EMAux_noaux$blimp_code)         # show generated BLIMP program code
blimp_print_psr(EMAux_noaux)        # PSR section only
blimp_print_focal(EMAux_noaux)      # focal outcome section only
rblimp::output(EMAux_noaux$fit)     # full output

#save focal outcome section
capture.output(
  blimp_print_focal(EMAux_noaux),
  file = file.path(output_path, "EMAux_noaux_focal_output.txt")
)

#save full output
capture.output(
  rblimp::output(EMAux_noaux$fit),
  file = file.path(output_path, "EMAux_noaux_full_output.txt")
)

#######################################################
### multilevel model with 13 auxiliary variables
#######################################################
EMAux_13aux <- EMAuxiliary(
  negaff ~  physact1*physact1.mean + (1 + physact1 | idn),
  data = data1,
  id = 'idn',
  aux = c("std_negaff", "lag_missing", "std_physact1", "lag_posaff", "age",
          "lag_physact2", "lag_negaff_ce", "ACC_NSteps",
          "lag_delay", "day_n", "BMI",
          "trigger_hr_sin", "weekend"),
  ordinal = c("weekend", "lag_missing"),
  center_group = 'physact1',
  center_grand = 'physact1.mean', 
  seed = 90291,
  burn = 20000,
  iter = 10000,
  chains = 12
)
cat(EMAux_13aux$blimp_code)         # show generated BLIMP program code
blimp_print_psr(EMAux_13aux)        # PSR section only
blimp_print_focal(EMAux_13aux)      # focal outcome section only
rblimp::output(EMAux_13aux$fit)     # full output

#save focal outcome section
capture.output(
  blimp_print_focal(EMAux_13aux),
  file = file.path(output_path, "EMAux_withaux_focal_output.txt")
)

#save full output
capture.output(
  rblimp::output(EMAux_13aux$fit),
  file = file.path(output_path, "EMAux_withaux_full_output.txt")
)


library(googledrive)
library(dplyr)

read_drive_csv <- function(folder_id, filename, sep) {
  files <- drive_ls(as_id(folder_id))
  target_file <- files %>% filter(name == filename)
  
  if (nrow(target_file) == 0) {
    stop(paste("File", filename, "not found in folder ID", folder_id))
  }
  
  temp_file <- tempfile(fileext = ".csv")
  drive_download(as_id(target_file$id), path = temp_file, overwrite = TRUE)
  data <- read.csv(temp_file, sep = sep, stringsAsFactors = FALSE)
  unlink(temp_file)
  
  return(data)
}

# STEP 1: Get IC folder in Shared Drive
ic_folder <- drive_ls(
  path = NULL,
  shared_drive = "Watershed_Biogeochemistry_Lab",
  type = "folder"
) %>% filter(name == "IC")

# STEP 2: Get IC_DATA
ic_data_folder <- drive_ls(as_id(ic_folder$id), type = "folder") %>%
  filter(name == "IC_DATA")

# STEP 3: Get project folder inside IC_DATA
project_path <- drive_ls(as_id(ic_data_folder$id), type = "folder") %>%
  filter(name == project_folder)

# STEP 4: Navigate to Standards_Runs > Anion and Cation folders
standards_runs_folder <- drive_ls(as_id(ic_data_folder$id), type = "folder") %>%
  filter(name == "Standards_Runs")

anion_cal_folder <- drive_ls(as_id(standards_runs_folder$id), type = "folder") %>%
  filter(name == "Anion")

cation_cal_folder <- drive_ls(as_id(standards_runs_folder$id), type = "folder") %>%
  filter(name == "Cation")


## if your personal email was linked to R, you need to delete that link and establish with the email that has access to Watershed Biogeochemistry Lab google drive 
# run this chuck of code first, then go back to main script to try again.
# 
# ###Disconnect links to reset
# drive_deauth()  # disconnect
# # Delete any cached token directories that might be interfering
# unlink(".secrets", recursive = TRUE)
# unlink("~/.R/gargle", recursive = TRUE)
# unlink("~/.cache/gargle", recursive = TRUE)
# 
# #Reset the link and make sure the Drive packages is associated with boisestate.edu account
# drive_auth(
#   cache = ".secrets",
#   email = "kyleformigli@u.boisestate.edu",
#   scopes = "https://www.googleapis.com/auth/drive"
# )
# 
# #checks to make sure user has access to the shared drive
# drive_ls(shared_drive = "Watershed_Biogeochemistry_Lab")
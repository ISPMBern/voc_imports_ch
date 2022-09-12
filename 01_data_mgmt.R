## Function loading of cleaned variables
# check initialize file if you installed all packages
path_name <- rstudioapi::getSourceEditorContext()$path
path_name <- gsub("01_data_mgmt.R","",path_name)
setwd(path_name)
source(paste0(path_name,"R/imports_000_initialize.R"))

# save data generated (is needed to generate rmd output file):
save.image(paste0("importVoC_project_",Sys.Date(),".RData"))

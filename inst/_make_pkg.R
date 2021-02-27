PKG <- "wbacon"
if (.Platform$OS.type == "unix") { 
	PKG_SOURCE <- "/mnt/c/my/code/cbacon"
	PKG_ROOT <- "/mnt/c/my/tmp"
} else {
	PKG_SOURCE <- "C:/My/code/cbacon"
	PKG_ROOT <- "C:/My/tmp"
}
setwd(PKG_ROOT)
#-------------------------------------------------------------------------------
# create folder 
if (dir.exists(PKG)) {
   unlink(PKG, recursive = TRUE)   
}

dir.create(PKG)
setwd(paste0(PKG_ROOT, "/", PKG))
folders <- c("man", "R", "src", "inst", "vignettes", "tests")
for (i in 1:length(folders)) {
   dir.create(folders[i])
}

#-------------------------------------------------------------------------------
# copy Rd-files
Rd_files <- list.files(paste0(PKG_SOURCE, "/man"), pattern = "\\.Rd$", 
   full.names = TRUE)
file.copy(Rd_files, paste0(PKG_ROOT, "/", PKG, "/man"), overwrite = TRUE)

# copy R-files
R_files <- list.files(paste0(PKG_SOURCE, "/R"), pattern = "\\.R$", 
   full.names = TRUE)
file.copy(R_files, paste0(PKG_ROOT, "/", PKG, "/R"), overwrite = TRUE)

# copy src-files
src_files <- list.files(paste0(PKG_SOURCE, "/src"), pattern = "\\.c$", 
   full.names = TRUE)
src_files <- c(src_files, list.files(paste0(PKG_SOURCE, "/src"), 
   pattern = "\\.h$", full.names = TRUE))
file.copy(src_files, paste0(PKG_ROOT, "/", PKG, "/src"), overwrite = TRUE)
#
file.copy(paste0(PKG_SOURCE, "/src/Makevars"), paste0(PKG_ROOT, "/", PKG, 
   "/src"), overwrite = TRUE)

# data-files 
# data_files <- list.files(paste0(PKG_SOURCE, "/data"), full.names = TRUE)
# file.copy(data_files, paste0(PKG_ROOT, "/", PKG, "/data"), overwrite = TRUE)

# inst-files 
inst_files <- list.files(paste0(PKG_SOURCE, "/inst"), full.names = TRUE)
file.copy(inst_files, paste0(PKG_ROOT, "/", PKG, "/inst"), overwrite = TRUE)

# vignettes-files 
vignettes_files <- list.files(paste0(PKG_SOURCE, "/vignettes"), 
   full.names = TRUE)
file.copy(vignettes_files, paste0(PKG_ROOT, "/", PKG, "/vignettes"), 
   overwrite = TRUE)

# tests-files 
tests_files <- list.files(paste0(PKG_SOURCE, "/tests"), 
   full.names = TRUE)
file.copy(tests_files, paste0(PKG_ROOT, "/", PKG, "/tests"), 
   overwrite = TRUE)

# NAMESPACE, DESCRIPTION and LICENSE
file.copy(paste0(PKG_SOURCE, "/", c("NAMESPACE", "DESCRIPTION", "LICENSE")), 
   paste0(PKG_ROOT, "/", PKG), overwrite = TRUE)

#-------------------------------------------------------------------------------
setwd(PKG_ROOT)
system("R CMD INSTALL --build wbacon")

# R CMD Rd2pdf wbacon # render pdf manual



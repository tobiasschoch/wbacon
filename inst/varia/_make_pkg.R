PKG <- "wbacon"
if (.Platform$OS.type == "unix") {
	PKG_SOURCE <- "/mnt/c/my/code"
	PKG_ROOT <- "/mnt/c/my/tmp"
} else {
	PKG_SOURCE <- "C:/My/code"
	PKG_ROOT <- "C:/My/tmp"
}
setwd(PKG_ROOT)

# delete package if it already exists
if (dir.exists(PKG))
    unlink(PKG, recursive = TRUE)

# copy entire directory (excl. files/folders with a leading dot, e.g. '.git')
dir.create(PKG)
pkg_files <- list.files(paste0(PKG_SOURCE, "/", PKG), full.names = TRUE)
file.copy(pkg_files, paste0(PKG_ROOT, "/", PKG), recursive = TRUE)

# clean src folder (remove binary files)
binary_files <- list.files(paste0(PKG_ROOT, "/", PKG, "/src"),
    pattern = "\\.o$|\\.dll$|\\.so$")
file.remove(paste0(PKG_ROOT, "/", PKG, "/src/", binary_files))

#-------------------------------------------------------------------------------
setwd(PKG_ROOT)

# fast install (only x64 arch; without html help files and vignette)
system("R CMD INSTALL wbacon --no-html --no-multiarch")

# install (without building the vignette)
#system("R CMD INSTALL wbacon")

# build the tar ball
#system("R CMD build wbacon")

# install from the tar ball
#system("R CMD INSTALL wbacon_0.2.tar.gz")

# install from the tar ball and generate *.zip-file
#system("R CMD INSTALL --build wbacon_0.2.tar.gz")

# render pdf manual
#system("R CMD Rd2pdf wbacon")

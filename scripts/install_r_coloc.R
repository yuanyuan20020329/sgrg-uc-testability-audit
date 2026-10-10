options(repos = c(CRAN = "https://cloud.r-project.org"))
lib <- Sys.getenv("R_LIB_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("libpaths:", .libPaths(), "\n")
for (p in c("coloc", "susieR")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, lib = lib, dependencies = FALSE)
  }
  cat(p, ":", as.character(packageVersion(p)), "\n")
}

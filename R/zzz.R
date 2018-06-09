.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) Rcpp::loadModule(m, what = TRUE)

  .jpackage(pkgname, lib.loc = libname)
  j = list.files(path= file.path(devtools::package_file(), 'inst','java'), pattern = ".jar$");
  if (length(j)==0) {print("Error: Java runtime library not found")
  } else {
    j=sort(j,decreasing=T); vm=file.path(devtools::package_file(), 'inst','java', j[1])
    .jinit (classpath=vm) # initialize java VM
  }

}

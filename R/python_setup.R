# python 'scipy' module I want to use in my package
ctef <- NULL

.onLoad <- function(libname, pkgname) {
  # delay load foo module (will only be loaded when accessed via $)
  # global reference to scipy (will be initialized in .onLoad)

    ctef <<- reticulate::import("ctef", delay_load = TRUE)

}

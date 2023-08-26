#'  global reference to ctef
#'
#'@description
#'`ctef` will be initialized in .onLoad.
#'
ctef <- NULL

#' Delay load ctef module
#'
#' @description
#' `.onLoad` delays loading ctef module (will only be loaded when accessed via $).
#'
#' @param libname Library name
#' @param pkgname Package name
.onLoad <- function(libname, pkgname) {
  # delay load foo module (will only be loaded when accessed via $)
  # global reference to scipy (will be initialized in .onLoad)

    ctef <<- reticulate::import("ctef", delay_load = TRUE)

}

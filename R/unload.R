.onUnload <- function (libpath) {
  library.dynam.unload("gpdIcm", libpath)
}

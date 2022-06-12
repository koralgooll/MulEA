.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to MulEA package!")
  assign("databaseConnection", NULL, envir = .GlobalEnv)
  assign("databaseLocalization", NULL, envir = .GlobalEnv)
}

#' @import ff methods parallel
#' @importFrom bit physical physical<-
NULL


.onAttach <- function(libname, pkgname) {
    packageStartupMessage("The symDMatrix package was supported by the National Institutes of Health (Grant: R01GM101219, R01GM099992).")
}

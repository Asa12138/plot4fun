#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom utils data head download.file untar combn
#' @importFrom stats runif rnorm
#' @importFrom grDevices colorRampPalette png pdf dev.off
#' @importFrom graphics plot.new text legend hist par
#' @importFrom pcutils lib_ps get_cols update_param download2
#' @import ggplot2
## usethis namespace: end
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

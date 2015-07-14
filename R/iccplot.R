#' Item Characteristic Curve
#'
#' The item characteristic curve is performed for the multidimensional polytomous
#' Rasch model or the continuous Rating Scale Model.
#'
#' The item characteristic curve (ICC) plots the response probability depending on person and item parameter.
#' For plotting the ICC, the object resulting from MPRM \code{\link{MPRM}} or CRSM \code{\link{CRSM}} is the input for the \code{iccplot} function.
#'
#' @aliases iccplot iccplot.CRSM iccplot.MPRM
#' @param object Object of class \code{CRSM} for ICC of the
#' CRSM or object of class \code{MPRM} for graphical model check of the MPRM
#' @param \dots \dots{}
#' @author Christine Hohensinn
#' @seealso \code{\link{LRT}} \code{\link{CRSM}}
#'
#' @rdname iccplot
#' @keywords item characteristic curve, item characteristic function
#' @examples
#'
#' #estimate CRSM for the first three items
#' data(analog)
#' res_cr <- CRSM(extraversion, low=-10, high=10)
#'
#' #ICC plot
#' iccplot(res_cr)
#'
#'
#' @export iccplot
iccplot <-
  function(object,...)UseMethod("iccplot")

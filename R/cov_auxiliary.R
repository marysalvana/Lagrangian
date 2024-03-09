#' Compute the univariate Lagrangian spatio-temporal stationary covariance function
#'
#' @description
#' \code{cov_uni_lagrangian} evaluates the univariate spatio-temporal stationary 
#' covariance function model based on the Lagrangian transport.
#'
#' @usage cov_uni_lagrangian(location, scale, vx, vy)
#'
#' @param location An \eqn{n \times 3} matrix of coordinates.
#' @param sigma2 A numeric constant indicating the spatial variance parameter.
#' @param scale A numeric constant indicating the spatial scale parameter.
#' @param nu A numeric constant indicating the spatial smoothness parameter.
#' @param vx A numeric constant indicating the x component of the transport velocity.
#' @param vy A numeric constant indicating the y component of the transport velocity.
#' @param nonstat 0, 1, 2 indicating the type of nonstationarity exhibited by the
#'        spatial field. Default is 0 (stationarity).
#'
#' @useDynLib Lagrangian, .registration=TRUE
#'
#' @return A matrix of dimension \eqn{n \times n}.
#'
#' @author Mary Lai Salvana \email{yourlainess@gmail.com}
#'
#' @examples
#'
#' library(dplyr)
#' library(fields)
#'
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#'
#' loc2d <- expand.grid(x, y) %>% as.matrix()
#' t <- seq(0, 10, by = 1)
#' loc3d <- cbind(rep(loc2d[, 1], each = length(t)), rep(loc2d[, 2], each = length(t)), t)
#'
#' cov_mat <- cov_uni_lagrangian(location = loc3d, sigma2 = 1, scale = 0.5, nu = 1, 
#'                               vx = 0.1001, vy = 0.1001)
#'
#'
#' @export
cov_uni_lagrangian <- function(location, sigma2, scale, nu, vx, vy, nonstat = 0){

  cov_mat <- CalculateCovUniLagrangian(x = location[, 1], y = location[, 2], t = location[, 3],
                                       sigma2 = sigma2, scale = scale, nu = nu,
                                       vx = vx, vy = vy, nonstat = nonstat)

  return(cov_mat)
}

CalculateCovUniLagrangian <- function (x, y, t, sigma2, scale, nu, vx, vy, nonstat) {

  n <- length(y)

  dist <- rep(0, n^2)

  res <- .C("CovUniLagrangian", x = as.double(x), y = as.double(y), t = as.double(t),
            sigma2 = as.double(sigma2), scale = as.double(scale), nu = as.double(nu),
            vx = as.double(vx), vx = as.double(vx),
            dist = as.double(dist), n = as.integer(n), nonstat = as.integer(nonstat))

  dist <- matrix(res$dist, n, n)
  return (dist)
}

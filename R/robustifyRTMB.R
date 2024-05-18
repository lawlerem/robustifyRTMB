#' Robustify likelihood components
#' 
#' Compute the robustified likelihood value
#' 
#' @param x A likelihood value
#' @param robustness A positive number. Larger numbers downweight components more.
#' @param robust_function Either "id" for the identity (non-robust), "ll" for log-logistic, or "ssh" for smooth semi-huber.
#' 
#' @return A numeric value giving the robustified likelihood component.
#' 
#' @export
robustify<- function(x, robustness, robust_function = "id") NULL

#' Compute the robust weight of a likelihood component
#' 
#' @return A numeric value giving the robust weight.
#' 
#' @describeIn robustify Compute the robust weight of a likelihood component
#' 
#' @export
robust_weight<- function(x, robustness, robust_function = "id") NULL

#' Compute the log-determinant of a matrix
#' 
#' @param m A square matrix
#' 
#' @return The log-determinant of m
#' 
#' @export
logdet<- function(m) NULL

.onLoad<- function(libname, pkgname) {
  .atomic_log_logistic<- RTMB::MakeTape(
    function(y) {
      x<- y[[1]]
      c<- y[[2]]
      ans<- log(1 + exp(x + c)) - log(1 + exp(c))
      return(ans)
    },
    c(0, 1)
  )$atomic()

  .orig_TapeConfig<- RTMB::TapeConfig()
  RTMB::TapeConfig(comparison = "tape")
  .atomic_smooth_semi_huber<- RTMB::MakeTape(
    function(y) {
      x<- y[[1]]
      c<- y[[2]]
      u<- (x + c) / c
      ans<- (x >= -c) * (x) + 
        (x < -c) * (c * log(u + sqrt(1 + u^2)) - c)
      return(ans)
    },
    c(0, 1)
  )$atomic()
  do.call(RTMB::TapeConfig, as.list(.orig_TapeConfig))

  robustify<- function(x, robustness, robust_function = "id") {
    c<- 1 / robustness
    ans<- switch(
      robust_function,
      "id" = x,
      "ll" = .atomic_log_logistic(c(x, c)),
      "ssh" = .atomic_smooth_semi_huber(c(x, c)),
      default = x
    )
    return(ans)
  }
  assign(
    "robustify",
    robustify,
    parent.env(environment())
  )

  robust_weight<- function(x, robustness, robust_function = "id") {
    c<- 1 / robustness
    ans<- switch(
      robust_function,
      "id" = 1,
      "ll" = .atomic_log_logistic$jacfun()(c(x, c))[1],
      "ssh" = .atomic_smooth_semi_huber$jacfun()(c(x, c))[1],
      default = 1
    )
    return(ans)
  }
  assign(
    "robust_weight",
    robust_weight,
    parent.env(environment())
  )

  logdet<- RTMB::ADjoint(
    function(x) {
        dim(x) <- rep(sqrt(length(x)), 2)
        determinant(x, log=TRUE)$modulus
    },
    function(x, y, dy) {
        dim(x) <- rep(sqrt(length(x)), 2)
        t(RTMB::solve(x)) * dy
    },
    name = "logdet"
  )
  assign(
    "logdet",
    logdet,
    parent.env(environment())
  )

  invisible()
}

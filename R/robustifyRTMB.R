#' Robustify likelihood components
#' 
#' Compute the robustified likelihood value
#' 
#' @param x A vector of log-likelihood values
#' @param robustness A positive number. Larger numbers downweight components more.
#' @param robust_function Either "id" for the identity (non-robust), "ll" for log-logistic, or "ssh" for smooth semi-huber.
#' 
#' @return The robustified log-likelihood values.
#' 
#' @export
robustify<- function(x, robustness, robust_function = "id") NULL

#' Compute the robust weight of a log-likelihood value
#' 
#' @return A vector of robust weights.
#' 
#' @describeIn robustify Compute the robust weight of log-likelihood values
#' 
#' @export
robust_weight<- function(x, robustness, robust_function = "id") NULL

.onLoad<- function(libname, pkgname) {
    .orig_TapeConfig<- RTMB::TapeConfig()
    RTMB::TapeConfig(comparison = "tape")

    .ll<- RTMB::MakeTape(
        function(y) {
            x<- y[[1]]
            c<- 1 / (y[[2]] + 0.002) # give upper bound to c to avoid overflow later
            ans<- log(1 + exp(x + c)) - log(1 + exp(c))
            ans<- (y[[2]] == 0) * x + (y[[2]] > 0) * ans
            return(ans)
        },
        c(0, 1)
    )$atomic()
    .ll_jac<- .ll$jacfun()

    .ssh<- RTMB::MakeTape(
        function(y) {
            x<- y[[1]]
            c<- 1 / (y[[2]] + .Machine$double.eps)
            u<- (x + c) / c
            ans<- (x >= -c) * (x) + 
                (x < -c) * (c * log(u + sqrt(1 + u^2)) - c)
            ans<- (y[[2]] == 0) * x + (y[[2]] > 0) * ans
            return(ans)
        },
        c(0, 1)
    )$atomic()
    .ssh_jac<- .ssh$jacfun()

    do.call(RTMB::TapeConfig, as.list(.orig_TapeConfig))

    robustify<- function(x, robustness, robust_function = "id") {
        ans<- RTMB::sapply(
            x,
            function(xx) {
                yy<- switch(
                    robust_function,
                    "id" = xx,
                    "ll" = .ll(c(xx, robustness)),
                    "ssh" = .ssh(c(xx, robustness)),
                    default = xx
                )
                return(yy)
            }
        )
        return(ans)
    }
    assign(
        "robustify",
        robustify,
        parent.env(environment())
    )

    robust_weight<- function(x, robustness, robust_function = "id") {
        ans<- RTMB::sapply(
            x,
            function(xx) {
                yy<- switch(
                    robust_function,
                    "id" = 1,
                    "ll" = .ll_jac(c(xx, robustness))[1],
                    "ssh" = .ssh_jac(c(xx, robustness))[1],
                    default = 1
                )
                return(yy)
            }
        )
        return(ans)
    }
    assign(
        "robust_weight",
        robust_weight,
        parent.env(environment())
    )

    invisible()
}

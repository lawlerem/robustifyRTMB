#' Fit a random effects model with a robustified likelihood
#' 
#' @param func 
#'     See RTMB::MakeADFun
#' @param parameters 
#'     See RTMB::MakeADFun
#' @param random 
#'     See RTMB::MakeADFun
#' @param smooth 
#'     A character naming the random effects that should be smoothed. 
#' @param nodes 
#'     A matrix giving the locations of the spline nodes
#' @param robust_schedule 
#'     The sequence of robustness values to iterate through
#' @param map 
#'     See RTMB::MakeADFun
#' @param bandwidth 
#'     The bandwidth to use when smoothing the spline between iterations
#' @param silent
#'     Should the optimization tracing be suppressed?
#' 
#' @return 
#'     A list containing:
#'   * par 
#'       The estimated parameter list
#'   * obj 
#'       The objective function from RTMB::MakeADFun
#'   * opt 
#'       The optimizer output from nlminb
#'   * sdr 
#'       The standard errors from RTMB::sdreport
#' 
#' @export
robustly_optimize<- function(
        func,
        parameters,
        random,
        smooth,
        nodes,
        robust_schedule = 0,
        map = list(),
        bandwidth = 1e-6,
        silent = TRUE
    ) {
    weight_fun<- (function(d, bandwidth) exp( -(d / bandwidth) )) |>
        Vectorize("d")
    d<- nodes |> stats::dist(diag = TRUE, upper = TRUE) |> as.matrix()
    within<- d <= (2 * bandwidth)
    weights<- 0 * d
    weights[within]<- d[within] |> weight_fun(bandwidth)
    weights<- weights |> 
        sweep(
            MARGIN = 1,
            STATS = weights |> rowSums(),
            FUN = "/"
        )
    tighten<- function(node_values) weights %*% cbind(node_values)

    if( "robustness" %in% names(parameters) ) {
        warning(
            paste0(
                "'robustness' is a reserved parameter for robustly_optimize. ",
                "Replacing it from the supplied parameter list."
            )
        )
        parameters$robustness<- NULL
    } else {}
    parameters<- c(
            parameters,
            list(
                robustness = robust_schedule[[1]]
            )
        ) |> 
        utils::as.relistable()
    if( "robustness" %in% names(map) ) map$robustness<- NULL
    robust_map<- list(robustness = as.factor(NA))
    spline_map<- list(
        node_values = NA |> 
            rep(parameters[[smooth]] |> length()) |>
            as.factor()
    )
    names(spline_map)<- smooth
    spline_map<- c(
        map[names(map) != smooth],
        spline_map
    )

    parameter_tracing<- robust_schedule |>
        lapply(
            function(r) {
                return(
                    list(
                        robustness = r,
                        parameters = parameters,
                        opt = NA
                    )
                )
            }
        )

    if( !silent & length(robust_schedule) > 0) cat("Creating taped function.\n")
    f<- func
    if( length(robust_schedule) > 1) f<- f |> RTMB::MakeTape(parameters)
    if( !silent ) cat(
            paste0(
                "(", parameters$robustness,
                " -> ", max(robust_schedule),
                ") Starting optimization."
            )
        )
    obj<- f |> 
        RTMB::MakeADFun(
            parameters,
            map = c(map, robust_map),
            random = random,
            silent = TRUE
        )
    opt<- suppressWarnings(with(obj, stats::nlminb(par, fn, gr)))
    fitpar<- obj$env$parList()
    parameter_tracing[[1]]$parameters<- fitpar
    parameter_tracing[[1]]$opt<- opt
    for( i in seq_along(robust_schedule)[-1] ) {
        r<- robust_schedule[i]
        fitpar$robustness<- r
        # 1.) Tighten spline estimate to remove pull from outliers
        if( !silent ) cat(
                paste0(
                    "\r\033[K(", fitpar$robustness,
                    " -> ", max(robust_schedule),
                    ") Tightening random effects. "
                )
            )
        fitpar[[smooth]]<- fitpar[[smooth]] |> tighten()

        # 2.) Re-estimate parameters holding spline fixed
        if( !silent ) cat(paste0("Re-estimating fixed parameters. "))
        obj<- f |> 
            RTMB::MakeADFun(
                fitpar,
                map = c(spline_map, robust_map),
                random = random,
                silent = TRUE
            )
        opt<- suppressWarnings(with(obj, stats::nlminb(par, fn, gr)))
        fitpar<- obj$env$parList()

        # 3.) Re-estimate spline and parameters
        if( !silent ) cat(paste0("Re-estimating all parameters."))
        if( !silent & (i == length(robust_schedule)) ) cat("\n")
        obj<- f |>
            RTMB::MakeADFun(
                fitpar,
                map = c(map, robust_map),
                random = random,
                silent = TRUE
            )
        opt<- suppressWarnings(with(obj, stats::nlminb(par, fn, gr)))
        fitpar<- obj$env$parList()

        parameter_tracing[[i]]$parameters<- fitpar
        parameter_tracing[[i]]$opt<- opt
    }

    convergences<- c |> 
        do.call(parameter_tracing |> lapply(function(x) x$opt$convergence))
    if( (convergences == 0) |> any() ) {
        last_converged<- (convergences == 0) |> which() |> max()
    } else {
        warning("Robust optimization did not converge.")
        last_converged<- convergences |> length()
    }
    if( !silent & (last_converged != length(convergences)) ) {
        cat(
            paste0(
                "Target robustness did not successfully converge. Using robustness ",
                parameter_tracing[[last_converged]]$robustness,
                " instead.\n"
            )
        )
    } else {}
    fitpar<- parameter_tracing[[last_converged]]$parameters
    opt<- parameter_tracing[[last_converged]]$opt

    if( !silent ) cat("Re-creating ADFun. ")
    if( length(robust_schedule) > 1 ) {
        obj<- func |>
            RTMB::MakeADFun(
                fitpar,
                map = c(map, robust_map),
                random = random,
                silent = TRUE
            )
    } else {}
    
    if( !silent ) cat("Computing sdreport.\n")
    sdr<- obj |> 
        RTMB::sdreport(
            opt$par,
            getJointPrecision = TRUE
        )
    
    return(
        list(
            par = fitpar,
            obj = obj,
            opt = opt,
            sdr = sdr
        )
    )
}





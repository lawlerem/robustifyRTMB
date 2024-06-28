#' Fit a random effects model with a robustified likelihood
#' 
#' @param func See RTMB::MakeADFun
#' @param parameters See RTMB::MakeADFun
#' @param random See RTMB::MakeADFun
#' @param smooth A character naming the random effects that should be smoothed
#' @param nodes A matrix giving the locations of the spline nodes
#' @param robust_schedule The sequence of robustness values to iterate through
#' @param map See RTMB::MakeADFun
#' @param bandwidth The bandwidth to use when smoothing the spline between iterations
#' @param inner.trace Should the inner tracing (nlminb) be printed?
#' @param outer.trace Should the outer tracing (iterative optimization) be printed?
#' 
#' @return A list containing the MakeADFun, the nlminb, and sdreport outputs
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
        inner.trace = FALSE,
        outer.trace = TRUE
    ) {
    weight_fun<- Vectorize(
        function(d, bandwidth) exp( -(d / bandwidth) ),
        "d"
    )
    d<- dist(
        nodes,
        diag = TRUE,
        upper = TRUE
    )
    d<- as.matrix(d)
    within<- d <= (2 * bandwidth)
    weights<- 0 * d
    weights[within]<- weight_fun(d[within], bandwidth)
    weights<- sweep(
        weights,
        MARGIN = 1,
        STATS = rowSums(weights),
        FUN = "/"
    )
    tighten<- function(node_values) weights %*% cbind(node_values)

    if( "robustness" %in% names(parameters) ) {
        warning("'robustness' is a reserved parameter for robustly_optimize. Replacing it from the supplied parameter list.")
        parameters$robustness<- NULL
    } else {}
    parameters<- c(
        parameters,
        list(
            robustness = robust_schedule[[1]]
        )
    )
    parameters<- as.relistable(parameters)
    if( "robustness" %in% names(map) ) map$robustness<- NULL
    robust_map<- list(
        robustness = as.factor(NA)
    )
    spline_map<- list(
        node_values = as.factor(
            rep(
                NA,
                length(parameters[[smooth]])
            )
        )
    )
    names(spline_map)<- smooth
    spline_map<- c(
        map[names(map) != smooth],
        spline_map
    )

    if( outer.trace & length(robust_schedule) > 0) cat("Creating taped function.\n")
    if( length(robust_schedule) == 1) {
        f<- func
    } else {
        f<- RTMB::MakeTape(
            func,
            parameters
        )
    }
    if( outer.trace ) cat(paste0("(", parameters$robustness, ") Starting optimization.\n"))
    obj<- RTMB::MakeADFun(
        f,
        parameters,
        map = c(map, robust_map),
        random = random,
        silent = !inner.trace
    )
    opt<- nlminb(
        obj$par,
        obj$fn,
        obj$gr
    )
    fitpar<- obj$env$parList()
    for( r in robust_schedule[-1] ) {
        fitpar$robustness<- r
        # 1.) Tighten spline estimate to remove pull from outliers
        if( outer.trace ) cat(paste0("(", fitpar$robustness, ") Tightening random effects. "))
        fitpar[[smooth]]<- tighten(fitpar[[smooth]])

        # 2.) Re-estimate parameters holding spline fixed
        if( outer.trace ) cat(paste0("Re-estimating fixed parameters. "))
        obj<- RTMB::MakeADFun(
            f,
            fitpar,
            map = c(spline_map, robust_map),
            random = random,
            silent = !inner.trace
        )
        opt<- nlminb(
            obj$par,
            obj$fn,
            obj$gr
        )
        fitpar<- obj$env$parList()

        # 3.) Re-estimate spline and parameters
        if( outer.trace ) cat(paste0("Re-estimating all parameters.\n"))
        obj<- RTMB::MakeADFun(
            f,
            fitpar,
            map = c(map, robust_map),
            random = random,
            silent = !inner.trace
        )
        opt<- nlminb(
            obj$par,
            obj$fn,
            obj$gr
        )
        fitpar<- obj$env$parList()
    }

    if( outer.trace ) cat("Re-creating ADFun. ")
    if( length(robust_schedule) > 1 ) {
        obj<- MakeADFun(
            func,
            fitpar,
            map = c(map, robust_map),
            random = random,
            silent = !inner.trace
        )
    } else {}
    
    if( outer.trace ) cat("Computing sdreport.\n")
    sdr<- sdreport(
        obj,
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





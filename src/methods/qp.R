# CIFresolve 0.1.1 quadratic-programming reconstruction.
CIFRESOLVE_IMPLEMENTATION <- "CIFresolve::KM_resolve"
CIFRESOLVE_VERSION <- as.character(utils::packageVersion("CIFresolve"))
CIFRESOLVE_GIT_COMMIT <- "2295c0d32b3bb5595fb38a32030f773daa47091d"
CIFRESOLVE_OPT_METHOD <- "approx"

.cifresolve_configuration <- function(configuration) {
  allowed <- c(
    "risk_ticks", "risk_ticks_deaths", "risk_only", "risk_deaths",
    "ticks_only", "ticks_deaths"
  )
  configuration <- match.arg(configuration, allowed)
  list(
    name = configuration,
    use_risk_table = startsWith(configuration, "risk_"),
    use_ticks = grepl("ticks", configuration, fixed = TRUE),
    use_deaths = grepl("deaths", configuration, fixed = TRUE)
  )
}

.prepare_cifresolve_input <- function(digi_dt, risk_table, arm, curve_id,
                                      censor_marks = NULL, total_events = NULL,
                                      configuration = "risk_ticks",
                                      maxdp = 8L) {
  cfg <- .cifresolve_configuration(configuration)

  xy <- data.table::as.data.table(.get_curve_xy(digi_dt, curve_id))
  data.table::setnames(xy, c("x", "y"), c("time", "surv"))
  xy <- xy[is.finite(time) & is.finite(surv) & time >= 0]
  xy[, `:=`(
    time = round(as.numeric(time), maxdp),
    surv = pmin(pmax(as.numeric(surv), 0), 1)
  )]
  xy <- xy[, .(surv = min(surv, na.rm = TRUE)), by = time]
  data.table::setorder(xy, time)
  xy[, surv := cummin(surv)]

  rt_full <- data.table::as.data.table(.get_arm_risk(risk_table, arm))
  data.table::setnames(rt_full, c("trisk", "nrisk"), c("time", "nrisk"))
  rt_full <- rt_full[is.finite(time) & !is.na(nrisk) & time >= 0]
  rt_full[, `:=`(time = round(as.numeric(time), maxdp), nrisk = as.integer(nrisk))]
  rt_full <- rt_full[, .(nrisk = min(nrisk, na.rm = TRUE)), by = time]
  data.table::setorder(rt_full, time)
  rt_full[, nrisk := cummin(nrisk)]

  if (!nrow(rt_full) || rt_full$time[1L] != 0) {
    stop("The number-at-risk table must start at time 0")
  }
  n0 <- as.integer(rt_full$nrisk[1L])
  tau <- max(rt_full$time)

  xy <- xy[time <= tau]
  if (!nrow(xy)) stop("No digitized KM coordinates remain in the reconstruction horizon")
  if (xy$time[1L] > 0) {
    xy <- data.table::rbindlist(list(data.table::data.table(time = 0, surv = 1), xy))
  }
  xy[time == 0, surv := 1]
  xy <- xy[, .(surv = min(surv, na.rm = TRUE)), by = time]
  data.table::setorder(xy, time)
  xy[, surv := cummin(surv)]

  if (cfg$use_risk_table) {
    rt <- data.table::copy(rt_full)
  } else {
    rt <- data.table::data.table(time = 0, nrisk = n0)
  }

  marks <- sort(unique(as.numeric(censor_marks %||% numeric(0))))
  marks <- marks[is.finite(marks) & marks >= 0 & marks <= tau]
  ticks <- if (cfg$use_ticks && length(marks)) marks else NULL

  ndeath <- NULL
  if (cfg$use_deaths) {
    ndeath <- .get_total_events(total_events, arm)
    if (is.null(ndeath) || length(ndeath) != 1L || !is.finite(ndeath)) {
      stop("A finite arm-specific total event count is required for ", configuration)
    }
    ndeath <- as.integer(ndeath)
  }

  list(
    S = list(time = as.numeric(xy$time), surv = as.numeric(xy$surv)),
    t.risk = as.numeric(rt$time),
    n.risk = as.integer(rt$nrisk),
    ticks = ticks,
    ndeath = ndeath,
    n0 = n0,
    tau = tau,
    n_curve_points = nrow(xy),
    n_risk_points = nrow(rt),
    n_unique_ticks = length(ticks %||% numeric(0)),
    configuration = configuration
  )
}

.validate_cifresolve_ipd <- function(ipd, n0) {
  if (!is.data.frame(ipd)) stop("CIFresolve::make_data() did not return a data frame")
  if (!nrow(ipd)) stop("CIFresolve::make_data() returned an empty reconstruction")
  if (any(!is.finite(ipd$time)) || any(ipd$time < 0)) {
    stop("CIFresolve returned invalid follow-up times")
  }
  if (any(!ipd$status %in% c(0L, 1L))) {
    stop("CIFresolve returned invalid event indicators")
  }
  invisible(TRUE)
}

.reconstruct_one_cifresolve <- function(digi_dt, risk_table, arm, curve_id,
                                        censor_marks = NULL, total_events = NULL,
                                        configuration = getOption(
                                          "km2ipd.cifresolve_configuration", "risk_ticks"
                                        ),
                                        strict_tick = getOption(
                                          "km2ipd.cifresolve_strict_tick", TRUE
                                        ),
                                        strict_dec = getOption(
                                          "km2ipd.cifresolve_strict_dec", TRUE
                                        ),
                                        constr_tol = getOption(
                                          "km2ipd.cifresolve_constr_tol", 1e-8
                                        ),
                                        cen_penalty = getOption(
                                          "km2ipd.cifresolve_cen_penalty", 1e-3
                                        ),
                                        maxdp = 8L) {
  prepared <- .prepare_cifresolve_input(
    digi_dt = digi_dt, risk_table = risk_table, arm = arm, curve_id = curve_id,
    censor_marks = censor_marks, total_events = total_events,
    configuration = configuration, maxdp = maxdp
  )

  captured_warnings <- character(0)
  fit <- withCallingHandlers(
    CIFresolve::KM_resolve(
      S = prepared$S,
      t.risk = prepared$t.risk,
      n.risk = prepared$n.risk,
      ndeath = prepared$ndeath,
      ticks = prepared$ticks,
      optmethod = CIFRESOLVE_OPT_METHOD,
      control = CIFresolve::control.CIFresolve(
        strict_tick = isTRUE(strict_tick),
        strict_dec = isTRUE(strict_dec),
        cen_penalty = as.numeric(cen_penalty),
        constr_tol = as.numeric(constr_tol),
        trace = 0,
        maxdp = as.integer(maxdp)
      )
    ),
    warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  raw_ipd <- CIFresolve::make_data(fit)
  ipd <- data.table::data.table(
    time = as.numeric(raw_ipd$time),
    status = as.integer(raw_ipd$event),
    arm = arm
  )
  .validate_cifresolve_ipd(ipd, prepared$n0)
  data.table::setorder(ipd, time, status)

  qp_object <- attr(fit, "object")
  diagnostics <- data.table::data.table(
    arm = arm,
    configuration = prepared$configuration,
    n_curve_points = prepared$n_curve_points,
    n_risk_points = prepared$n_risk_points,
    n_unique_ticks = prepared$n_unique_ticks,
    supplied_total_events = prepared$ndeath %||% NA_integer_,
    expected_sample_size = prepared$n0,
    reconstructed_sample_size = nrow(ipd),
    sample_size_difference = nrow(ipd) - prepared$n0,
    reconstructed_events = sum(ipd$status == 1L),
    reconstructed_censors = sum(ipd$status == 0L),
    continuous_objective = if (is.list(qp_object) && !is.null(qp_object$value)) {
      as.numeric(qp_object$value)
    } else NA_real_,
    integer_objective = as.numeric(attr(fit, "int_obj") %||% NA_real_),
    warning = if (length(captured_warnings)) {
      paste(unique(captured_warnings), collapse = " | ")
    } else NA_character_
  )

  list(ipd = ipd, diagnostics = diagnostics, fit = fit)
}

reconstruct_twoarm_cifresolve <- function(digi, risk_table,
                                          curve_map = c(Control = 1, Treatment = 2),
                                          total_events = NULL,
                                          censor_marks_by_arm = NULL,
                                          configuration = getOption(
                                            "km2ipd.cifresolve_configuration", "risk_ticks"
                                          )) {
  digi_dt <- if (is.list(digi) && !is.null(digi$data)) digi$data else digi
  digi_dt <- data.table::as.data.table(digi_dt)
  arms <- c("Control", "Treatment")
  fits <- lapply(arms, function(arm) {
    marks <- if (!is.null(censor_marks_by_arm)) {
      censor_marks_by_arm[[arm]] %||% numeric(0)
    } else numeric(0)
    .reconstruct_one_cifresolve(
      digi_dt = digi_dt,
      risk_table = risk_table,
      arm = arm,
      curve_id = unname(curve_map[[arm]]),
      censor_marks = marks,
      total_events = total_events,
      configuration = configuration
    )
  })
  ipd <- data.table::rbindlist(lapply(fits, `[[`, "ipd"), use.names = TRUE)
  ipd <- .fix_arm(ipd)
  list(
    method = "QP",
    ipd = ipd,
    ipd_control = fits[[1L]]$ipd,
    ipd_treat = fits[[2L]]$ipd,
    diagnostics = data.table::rbindlist(lapply(fits, `[[`, "diagnostics"), fill = TRUE)
  )
}

.reconstruct_one_qp <- function(digi_dt, risk_table, arm, curve_id,
                                 censor_marks = NULL, drop_tol = 1e-8,
                                 risk_align = TRUE, eps_time = 1e-7) {
  ans <- .reconstruct_one_cifresolve(
    digi_dt = digi_dt,
    risk_table = risk_table,
    arm = arm,
    curve_id = curve_id,
    censor_marks = censor_marks,
    total_events = getOption("km2ipd.cifresolve_total_events", NULL),
    configuration = getOption("km2ipd.cifresolve_configuration", "risk_ticks")
  )
  ans$ipd
}

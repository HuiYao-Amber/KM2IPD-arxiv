# Survival simulation, numeric KM digitization, and reconstruction metrics.
.require_pkgs <- function(pkgs, stop_if_missing = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L && stop_if_missing) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      "\nPlease install with install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), "))",
      call. = FALSE
    )
  }
  invisible(length(missing) == 0L)
}

.require_pkgs(c("survival", "data.table"), stop_if_missing = TRUE)
suppressPackageStartupMessages({
  library(survival)
  library(data.table)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

.safe_num <- function(x) suppressWarnings(as.numeric(x))

.as_dt <- function(x) data.table::as.data.table(x)

.fix_arm <- function(ipd) {
  ipd <- .as_dt(ipd)
  stopifnot(all(c("time", "status", "arm") %in% names(ipd)))
  ipd[, arm := as.character(arm)]
  ipd[arm %in% c("0", "control", "Control", "C", "ctrl"), arm := "Control"]
  ipd[arm %in% c("1", "treat", "Treatment", "T", "trt"), arm := "Treatment"]
  ipd[, arm := factor(arm, levels = c("Control", "Treatment"))]
  ipd[, time := as.numeric(time)]
  ipd[, status := as.integer(status)]
  ipd[]
}

.weibull_scale_from_mean <- function(mean_time, shape) {
  mean_time / gamma(1 + 1 / shape)
}

make_weibull_ph_scenario <- function(id, shape, hr, mean_time = 10) {
  scale_ctrl <- .weibull_scale_from_mean(mean_time, shape)
  scale_trt  <- scale_ctrl / (hr^(1 / shape))
  list(
    id = id,
    phase = "phase1",
    type = "weibull_ph",
    label = sprintf("Weibull PH: shape=%.2f, HR=%.2f", shape, hr),
    shape = shape,
    hr = hr,
    mean_time = mean_time,
    scale_ctrl = scale_ctrl,
    scale_trt = scale_trt
  )
}

make_phase1_scenario_library <- function(mean_time = 10) {
  list(
    W_shape1_HR067  = make_weibull_ph_scenario("W_shape1_HR067",  shape = 1.0, hr = 0.67, mean_time = mean_time),
    W_shape1_HR085  = make_weibull_ph_scenario("W_shape1_HR085",  shape = 1.0, hr = 0.85, mean_time = mean_time),
    W_shape1_HR100  = make_weibull_ph_scenario("W_shape1_HR100",  shape = 1.0, hr = 1.00, mean_time = mean_time),
    W_shape0p6_HR067 = make_weibull_ph_scenario("W_shape0p6_HR067", shape = 0.6, hr = 0.67, mean_time = mean_time),
    W_shape2_HR067  = make_weibull_ph_scenario("W_shape2_HR067",  shape = 2.0, hr = 0.67, mean_time = mean_time)
  )
}

make_piecewise_scenario <- function(id, label, control_haz, treat_haz, breaks = c(0, 6, 12), nominal_hr = NA_real_) {
  stopifnot(length(control_haz) == length(treat_haz), length(control_haz) == length(breaks))
  list(
    id = id,
    phase = "phase2",
    type = "piecewise_exp",
    label = label,
    breaks = breaks,
    control_haz = control_haz,
    treat_haz = treat_haz,
    nominal_hr = nominal_hr
  )
}

make_phase2_scenario_library <- function() {
  breaks <- c(0, 6, 12)
  list(
    PH_H1_HR067 = make_piecewise_scenario(
      id = "PH_H1_HR067",
      label = "Piecewise PH: constant hazard, HR=0.67",
      control_haz = c(0.080, 0.080, 0.080),
      treat_haz   = c(0.080, 0.080, 0.080) * 0.67,
      breaks = breaks,
      nominal_hr = 0.67
    ),
    PH_H2_HR085 = make_piecewise_scenario(
      id = "PH_H2_HR085",
      label = "Piecewise PH: increasing baseline hazard, HR=0.85",
      control_haz = c(0.045, 0.080, 0.120),
      treat_haz   = c(0.045, 0.080, 0.120) * 0.85,
      breaks = breaks,
      nominal_hr = 0.85
    ),
    PH_H2_HR100 = make_piecewise_scenario(
      id = "PH_H2_HR100",
      label = "Piecewise PH: increasing baseline hazard, HR=1.00",
      control_haz = c(0.045, 0.080, 0.120),
      treat_haz   = c(0.045, 0.080, 0.120),
      breaks = breaks,
      nominal_hr = 1.00
    ),
    NPHdelayed_H1_D1 = make_piecewise_scenario(
      id = "NPHdelayed_H1_D1",
      label = "NPH delayed benefit: early neutral/worse, late benefit",
      control_haz = c(0.080, 0.080, 0.080),
      treat_haz   = c(0.095, 0.075, 0.045),
      breaks = breaks,
      nominal_hr = NA_real_
    ),
    NPHdelayed_H2_D2 = make_piecewise_scenario(
      id = "NPHdelayed_H2_D2",
      label = "NPH delayed benefit: increasing baseline hazard",
      control_haz = c(0.045, 0.080, 0.120),
      treat_haz   = c(0.055, 0.075, 0.080),
      breaks = breaks,
      nominal_hr = NA_real_
    ),
    NPHcross_H1_C1 = make_piecewise_scenario(
      id = "NPHcross_H1_C1",
      label = "NPH crossing: early harm, late benefit",
      control_haz = c(0.060, 0.080, 0.100),
      treat_haz   = c(0.130, 0.075, 0.035),
      breaks = breaks,
      nominal_hr = NA_real_
    )
  )
}

get_scenario_library <- function() {
  c(make_phase1_scenario_library(), make_phase2_scenario_library())
}

.r_piecewise_exp <- function(n, hazards, breaks) {
  stopifnot(length(hazards) == length(breaks), all(hazards > 0))
  n <- as.integer(n)
  out <- numeric(n)
  for (i in seq_len(n)) {
    e <- stats::rexp(1, rate = 1)
    t <- breaks[1]
    for (j in seq_along(hazards)) {
      h <- hazards[j]
      start <- breaks[j]
      end <- if (j < length(breaks)) breaks[j + 1L] else Inf
      width <- end - start
      segment_H <- h * width
      if (e <= segment_H || !is.finite(width)) {
        out[i] <- start + e / h
        break
      } else {
        e <- e - segment_H
      }
    }
  }
  out
}

.generate_latent_event_times <- function(scn, n_control, n_treatment) {
  if (identical(scn$type, "weibull_ph")) {
    T_control <- stats::rweibull(as.integer(n_control), shape = scn$shape, scale = scn$scale_ctrl)
    T_treat   <- stats::rweibull(as.integer(n_treatment), shape = scn$shape, scale = scn$scale_trt)
  } else if (identical(scn$type, "piecewise_exp")) {
    T_control <- .r_piecewise_exp(as.integer(n_control), hazards = scn$control_haz, breaks = scn$breaks)
    T_treat   <- .r_piecewise_exp(as.integer(n_treatment), hazards = scn$treat_haz, breaks = scn$breaks)
  } else {
    stop("Unknown scenario type: ", scn$type)
  }
  list(T_control = T_control, T_treat = T_treat)
}

.make_uniforms <- function(n_control, n_treatment) {
  list(
    U_c_control = stats::runif(as.integer(n_control)),
    U_c_treat   = stats::runif(as.integer(n_treatment)),
    U_info_control = stats::runif(as.integer(n_control)),
    U_info_treat   = stats::runif(as.integer(n_treatment))
  )
}

.make_censor_times <- function(scale_c, censoring, target_censoring, T_control, T_treat, U) {
  censoring <- match.arg(censoring, c("random", "exp", "front", "back", "informative"))
  scale_c <- max(as.numeric(scale_c), 1e-8)

  if (censoring == "random") {
    C_control <- U$U_c_control * scale_c
    C_treat   <- U$U_c_treat * scale_c
  } else if (censoring == "front") {
    C_control <- stats::qbeta(U$U_c_control, shape1 = 0.7, shape2 = 2.0) * scale_c
    C_treat   <- stats::qbeta(U$U_c_treat,   shape1 = 0.7, shape2 = 2.0) * scale_c
  } else if (censoring == "back") {
    C_control <- stats::qbeta(U$U_c_control, shape1 = 2.0, shape2 = 0.7) * scale_c
    C_treat   <- stats::qbeta(U$U_c_treat,   shape1 = 2.0, shape2 = 0.7) * scale_c
  } else if (censoring == "exp") {
    lambda <- 2 / scale_c
    qtrunc <- function(u) stats::qexp(u * (1 - exp(-lambda * scale_c)), rate = lambda)
    C_control <- qtrunc(U$U_c_control)
    C_treat   <- qtrunc(U$U_c_treat)
  } else {
    p_inf <- min(max(target_censoring, 0.01), 0.85)
    C_control <- rep(Inf, length(T_control))
    C_treat   <- rep(Inf, length(T_treat))
    idx_c <- U$U_info_control < p_inf
    idx_t <- U$U_info_treat < p_inf
    C_control[idx_c] <- U$U_c_control[idx_c] * pmin(T_control[idx_c], scale_c)
    C_treat[idx_t]   <- U$U_c_treat[idx_t]   * pmin(T_treat[idx_t], scale_c)
  }
  list(C_control = C_control, C_treat = C_treat)
}

.realize_ipd_given_scale <- function(scale_c, censoring, target_censoring, T_control, T_treat, U, admin_time = NULL) {
  C <- .make_censor_times(scale_c, censoring, target_censoring, T_control, T_treat, U)
  admin_time <- admin_time %||% scale_c

  X_control <- pmin(T_control, C$C_control, admin_time)
  D_control <- as.integer(T_control <= C$C_control & T_control <= admin_time)

  X_treat <- pmin(T_treat, C$C_treat, admin_time)
  D_treat <- as.integer(T_treat <= C$C_treat & T_treat <= admin_time)

  ipd <- data.table::data.table(
    time = c(X_control, X_treat),
    status = c(D_control, D_treat),
    arm = c(rep("Control", length(X_control)), rep("Treatment", length(X_treat)))
  )
  ipd[, arm := factor(arm, levels = c("Control", "Treatment"))]
  data.table::setorder(ipd, arm, time)

  list(
    ipd = ipd,
    censor_prop = mean(ipd$status == 0),
    total_events = c(Control = sum(D_control == 1L), Treatment = sum(D_treat == 1L)),
    censor_times = C
  )
}

.calibrate_censor_scale <- function(target_censoring, censoring, T_control, T_treat, U,
                                    admin_time = NULL, tol = 1e-4, scale_hi_init = NULL,
                                    scale_cap = 1e5) {
  target_censoring <- max(min(target_censoring, 0.95), 0)

  f_root <- function(scale_c) {
    .realize_ipd_given_scale(scale_c, censoring, target_censoring, T_control, T_treat, U, admin_time)$censor_prop - target_censoring
  }

  lo <- 1e-6
  hi <- scale_hi_init %||% max(5, stats::quantile(c(T_control, T_treat), 0.90, na.rm = TRUE))

  flo <- suppressWarnings(f_root(lo))
  fhi <- suppressWarnings(f_root(hi))

  iter <- 0L
  while (is.finite(flo) && is.finite(fhi) && flo * fhi > 0 && iter < 60L && hi < scale_cap) {
    hi <- min(hi * 1.6, scale_cap)
    fhi <- suppressWarnings(f_root(hi))
    iter <- iter + 1L
  }

  if (is.finite(flo) && is.finite(fhi) && flo * fhi <= 0) {
    return(stats::uniroot(f_root, interval = c(lo, hi), tol = tol)$root)
  }

  grid <- exp(seq(log(1e-4), log(scale_cap), length.out = 100))
  vals <- vapply(grid, function(z) abs(f_root(z)), numeric(1))
  grid[which.min(vals)]
}

compute_pretty_axis <- function(ipd, multiple = 8L, fixed_x_end = NULL) {
  ipd <- .fix_arm(ipd)
  x_end_raw <- max(ipd$time, na.rm = TRUE)
  x_end_plot <- fixed_x_end %||% {
    x1 <- ceiling(as.numeric(x_end_raw)[1])
    as.numeric(ceiling(x1 / multiple) * multiple)
  }
  x_increment <- x_end_plot / multiple
  ticks <- seq(0, x_end_plot, by = x_increment)
  list(x_end_raw = x_end_raw, x_end_plot = x_end_plot, x_increment = x_increment, ticks = ticks)
}

make_risk_table <- function(ipd, times) {
  ipd <- .fix_arm(ipd)
  nrisk <- function(tt, tvec) sum(tvec >= tt)
  time_ctrl <- ipd[arm == "Control", time]
  time_trt <- ipd[arm == "Treatment", time]
  rt <- data.table::data.table(time = as.numeric(times))
  rt[, nrisk_control := vapply(time, nrisk, numeric(1), tvec = time_ctrl)]
  rt[, nrisk_treat := vapply(time, nrisk, numeric(1), tvec = time_trt)]
  rt[, nrisk_all := nrisk_control + nrisk_treat]
  rt[]
}

extract_censor_marks_by_arm <- function(ipd, round_digits = NULL) {
  ipd <- .fix_arm(ipd)
  marks <- split(ipd[status == 0L, time], ipd[status == 0L, arm])
  for (a in c("Control", "Treatment")) {
    if (is.null(marks[[a]])) marks[[a]] <- numeric(0)
    marks[[a]] <- sort(as.numeric(marks[[a]]))
    if (!is.null(round_digits)) marks[[a]] <- sort(unique(round(marks[[a]], round_digits)))
  }
  marks
}

simulate_single_trial <- function(scenario_id, scenario_lib = get_scenario_library(),
                                  n_control = 250, n_treatment = 250,
                                  censoring = c("random", "exp", "front", "back", "informative"),
                                  target_censoring = 0.30,
                                  admin_time = NULL,
                                  fixed_x_end = NULL,
                                  seed = NULL) {
  censoring <- match.arg(censoring)
  if (!is.null(seed)) set.seed(seed)
  scn <- scenario_lib[[scenario_id]]
  if (is.null(scn)) stop("Unknown scenario_id: ", scenario_id)

  latent <- .generate_latent_event_times(scn, n_control, n_treatment)
  U <- .make_uniforms(n_control, n_treatment)

  if (is.null(admin_time) && identical(scn$phase, "phase2")) admin_time <- 24
  scale_c <- .calibrate_censor_scale(
    target_censoring = target_censoring,
    censoring = censoring,
    T_control = latent$T_control,
    T_treat = latent$T_treat,
    U = U,
    admin_time = admin_time
  )
  realized <- .realize_ipd_given_scale(scale_c, censoring, target_censoring,
                                       latent$T_control, latent$T_treat, U, admin_time)
  ipd <- realized$ipd

  axis <- compute_pretty_axis(ipd, multiple = 8L, fixed_x_end = fixed_x_end %||% admin_time)
  risk_table <- make_risk_table(ipd, axis$ticks)

  list(
    scenario = scn,
    scenario_id = scenario_id,
    phase = scn$phase,
    ipd = ipd,
    risk_table = risk_table,
    axis = axis,
    censor_marks_by_arm = extract_censor_marks_by_arm(ipd),
    total_events = realized$total_events,
    censor_prop = realized$censor_prop,
    censor_scale = scale_c,
    admin_time = admin_time %||% scale_c,
    latent = latent
  )
}

.km_surv_at_times <- function(ipd_arm, times) {
  ipd_arm <- data.table::as.data.table(ipd_arm)
  if (nrow(ipd_arm) == 0L) return(rep(NA_real_, length(times)))
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd_arm)
  if (length(fit$time) == 0L) return(rep(1, length(times)))
  out <- summary(fit, times = times, extend = TRUE)$surv
  if (anyNA(out)) {
    for (i in seq_along(out)) {
      if (is.na(out[i])) out[i] <- if (i == 1L) 1 else out[i - 1L]
    }
  }
  out
}

digitize_km_numeric <- function(ipd_true, risk_table = NULL,
                                mode = c("event_times", "grid"),
                                n_grid = 200,
                                sd_S = 0.01,
                                include_risk_times = TRUE,
                                x_end_plot = NULL,
                                seed = NULL) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  ipd_true <- .fix_arm(ipd_true)
  axis <- compute_pretty_axis(ipd_true, fixed_x_end = x_end_plot)
  x_end <- x_end_plot %||% axis$x_end_plot

  make_one <- function(arm_label, curve_id) {
    dat <- ipd_true[arm == arm_label]
    if (mode == "event_times") {
      ev_times <- sort(unique(dat[status == 1L & time <= x_end, time]))
      times <- unique(c(0, ev_times, x_end))
      if (include_risk_times && !is.null(risk_table)) times <- unique(c(times, risk_table$time))
      times <- sort(times[is.finite(times) & times >= 0 & times <= x_end])
    } else {
      times <- seq(0, x_end, length.out = as.integer(n_grid))
      if (include_risk_times && !is.null(risk_table)) times <- sort(unique(c(times, risk_table$time)))
    }
    St <- .km_surv_at_times(dat, times)
    St_noisy <- St + stats::rnorm(length(St), mean = 0, sd = sd_S)
    St_noisy[times == min(times)] <- 1
    St_noisy <- pmin(pmax(St_noisy, 0), 1)
    St_noisy <- cummin(St_noisy)
    St_noisy[times == min(times)] <- 1
    data.table::data.table(time = times, St = St_noisy, curve = as.integer(curve_id))
  }

  dC <- make_one("Control", 1L)
  dT <- make_one("Treatment", 2L)
  digi <- data.table::rbindlist(list(dC, dT), use.names = TRUE)
  data.table::setorder(digi, curve, time)
  list(data = digi, axis = axis, mode = mode, sd_S = sd_S, x_end_plot = x_end)
}

.get_curve_xy <- function(digi_dt, curve_id) {
  dt <- .as_dt(digi_dt)
  d <- dt[curve == curve_id, .(time, St)]
  data.table::setorder(d, time)
  d[, St := pmin(pmax(St, 0), 1)]
  d[, St := cummin(St)]
  data.frame(x = d$time, y = d$St)
}

.get_arm_risk <- function(risk_table, arm) {
  rt <- .as_dt(risk_table)
  data.table::setorder(rt, time)
  nrisk <- if (arm == "Control") rt$nrisk_control else rt$nrisk_treat
  data.frame(trisk = rt$time, nrisk = as.integer(nrisk))
}

.get_total_events <- function(total_events, arm) {
  if (is.null(total_events)) return(NULL)
  if (is.list(total_events) && !is.null(total_events[[arm]])) return(as.integer(total_events[[arm]]))
  if (!is.null(names(total_events)) && arm %in% names(total_events)) return(as.integer(total_events[[arm]]))
  NULL
}

.extract_drop_bottoms <- function(digi_dt, curve_id, drop_tol = 1e-8) {
  dt <- .as_dt(digi_dt)
  d <- dt[curve == curve_id, .(time, St)]
  data.table::setorder(d, time)
  d <- d[, .(St = min(St, na.rm = TRUE)), by = time]
  data.table::setorder(d, time)
  d[, St := pmin(pmax(St, 0), 1)]
  d[, St := cummin(St)]
  idx <- which(diff(d$St) < -abs(drop_tol)) + 1L
  data.frame(t = d$time[idx], S = d$St[idx])
}

.ensure_ipdfromkm_namespace <- function() {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package dplyr is required by IPDfromKM::preprocess(). Please install.packages('dplyr').")
  }
  suppressPackageStartupMessages(library(dplyr))
  ns <- tryCatch(asNamespace("IPDfromKM"), error = function(e) NULL)
  if (!is.null(ns)) {
    try({
      if (!exists("n", envir = ns, inherits = FALSE)) {
        assign("n", dplyr::n, envir = ns)
      }
    }, silent = TRUE)
    try(utils::assignInNamespace("n", dplyr::n, ns = "IPDfromKM"), silent = TRUE)
  }
  if (!exists("n", envir = .GlobalEnv, inherits = FALSE)) {
    try(assign("n", dplyr::n, envir = .GlobalEnv), silent = TRUE)
  }
  invisible(TRUE)
}

.clean_xy_for_reconstruction <- function(xy) {
  xy <- data.table::as.data.table(xy)
  if ("x" %in% names(xy)) data.table::setnames(xy, "x", "time")
  if ("y" %in% names(xy)) data.table::setnames(xy, "y", "St")
  xy <- xy[is.finite(time) & is.finite(St), .(time = as.numeric(time), St = as.numeric(St))]
  xy[, St := pmin(pmax(St, 0), 1)]
  xy <- xy[, .(St = min(St, na.rm = TRUE)), by = time]
  data.table::setorder(xy, time)
  xy[, St := cummin(St)]
  if (nrow(xy) == 0L || xy$time[1] > 0 || xy$St[1] < 1) {
    xy <- data.table::rbindlist(list(data.table::data.table(time = 0, St = 1), xy), use.names = TRUE)
    xy <- xy[, .(St = min(St, na.rm = TRUE)), by = time]
    data.table::setorder(xy, time)
    xy[, St := cummin(St)]
  }
  as.data.frame(xy[, .(x = time, y = St)])
}

.make_ipdfromkm_step_xy <- function(xy) {
  xy <- data.table::as.data.table(xy)
  data.table::setnames(xy, c("x", "y"), c("time", "surv"))
  xy <- xy[is.finite(time) & is.finite(surv)]
  data.table::setorder(xy, time)
  xy[, surv := pmin(pmax(surv, 0), 1)]
  xy[, surv := cummin(surv)]
  if (nrow(xy) == 0L || xy$time[1] > 0 || xy$surv[1] < 1) {
    xy <- data.table::rbindlist(list(data.table::data.table(time = 0, surv = 1), xy), use.names = TRUE)
  }
  data.table::setorder(xy, time, -surv)
  xy <- unique(xy, by = c("time", "surv"))
  out <- data.table::data.table(time = numeric(0), surv = numeric(0))
  prev <- 1
  for (i in seq_len(nrow(xy))) {
    ti <- as.numeric(xy$time[i])
    si <- as.numeric(xy$surv[i])
    if (ti == 0) {
      out <- data.table::rbindlist(list(out, data.table::data.table(time = 0, surv = 1)), use.names = TRUE)
      prev <- min(prev, si)
      next
    }
    out <- data.table::rbindlist(list(
      out,
      data.table::data.table(time = ti, surv = prev),
      data.table::data.table(time = ti, surv = si)
    ), use.names = TRUE)
    prev <- min(prev, si)
  }
  out <- unique(out, by = c("time", "surv"))
  data.table::setorder(out, time, -surv)
  as.data.frame(out[, .(x = time, y = surv)])
}

.add_risk_time_anchors <- function(xy, rt) {
  xdt <- data.table::as.data.table(xy)
  data.table::setnames(xdt, c("x", "y"), c("time", "surv"))
  xdt <- xdt[is.finite(time) & is.finite(surv)]
  data.table::setorder(xdt, time, -surv)
  rtimes <- sort(unique(as.numeric(rt$trisk)))
  rtimes <- rtimes[is.finite(rtimes)]
  if (length(rtimes) == 0L) return(as.data.frame(xdt[, .(x = time, y = surv)]))
  max_x <- max(xdt$time, na.rm = TRUE)
  anchors <- data.table::data.table(time = numeric(0), surv = numeric(0))
  for (tt in rtimes) {
    if (!is.finite(tt) || tt < 0 || tt > max_x) next
    idx <- which(xdt$time <= tt)
    ss <- if (length(idx) == 0L) 1 else min(xdt$surv[idx], na.rm = TRUE)
    anchors <- data.table::rbindlist(list(anchors, data.table::data.table(time = tt, surv = ss)), use.names = TRUE)
  }
  out <- data.table::rbindlist(list(xdt, anchors), use.names = TRUE, fill = TRUE)
  out <- unique(out, by = c("time", "surv"))
  data.table::setorder(out, time, -surv)
  as.data.frame(out[, .(x = time, y = surv)])
}

.call_ipdfromkm_getIPD <- function(prep, armID, tot_ev = NULL) {
  fml <- names(formals(IPDfromKM::getIPD))
  arm_arg <- if ("armID" %in% fml) "armID" else if ("armind" %in% fml) "armind" else NULL
  if (is.null(arm_arg)) stop("Cannot identify IPDfromKM::getIPD arm argument name")

  make_call <- function(use_tot = TRUE) {
    args <- list(prep = prep)
    args[[arm_arg]] <- armID
    if (isTRUE(use_tot) && !is.null(tot_ev) && is.finite(tot_ev)) args[["tot.events"]] <- as.integer(tot_ev)
    do.call(IPDfromKM::getIPD, args)
  }

  errors <- character(0)
  out <- tryCatch(make_call(TRUE), error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL })
  if (!is.null(out)) return(out)
  out <- tryCatch(make_call(FALSE), error = function(e) { errors <<- c(errors, conditionMessage(e)); NULL })
  if (!is.null(out)) return(out)
  stop(paste(unique(errors), collapse = " | "))
}

.reconstruct_one_ipdfromkm <- function(digi_dt, risk_table, arm, curve_id, total_events = NULL) {
  if (!requireNamespace("IPDfromKM", quietly = TRUE)) stop("Package IPDfromKM is not installed")
  .ensure_ipdfromkm_namespace()
  xy0 <- .get_curve_xy(digi_dt, curve_id)
  xy_post <- .clean_xy_for_reconstruction(xy0)
  rt <- .get_arm_risk(risk_table, arm)
  data.table::setDT(rt)
  data.table::setorder(rt, trisk)
  n0 <- as.integer(rt$nrisk[which.min(rt$trisk)])
  armID <- if (arm == "Control") 0 else 1
  tot_ev <- .get_total_events(total_events, arm)

  xy_step <- .make_ipdfromkm_step_xy(xy_post)
  xy_candidates <- list(
    step_anchor = .add_risk_time_anchors(xy_step, rt),
    step = xy_step,
    post_anchor = .add_risk_time_anchors(xy_post, rt),
    post = xy_post
  )

  errors <- character(0)
  for (xy_name in names(xy_candidates)) {
    xy <- xy_candidates[[xy_name]]
    if (NROW(xy) < 5L) next
    prep_candidates <- list()

    prep1 <- tryCatch(
      IPDfromKM::preprocess(dat = xy, trisk = rt$trisk, nrisk = rt$nrisk, totalpts = n0, maxy = 1),
      error = function(e) { errors <<- c(errors, paste0(xy_name, "/risk-prep: ", conditionMessage(e))); NULL }
    )
    if (!is.null(prep1)) prep_candidates[[paste0(xy_name, "/risk")]] <- prep1

    prep2 <- tryCatch(
      IPDfromKM::preprocess(dat = xy, totalpts = n0, maxy = 1),
      error = function(e) { errors <<- c(errors, paste0(xy_name, "/nonly-prep: ", conditionMessage(e))); NULL }
    )
    if (!is.null(prep2)) prep_candidates[[paste0(xy_name, "/nonly")]] <- prep2

    for (prep_name in names(prep_candidates)) {
      out <- tryCatch(
        .call_ipdfromkm_getIPD(prep = prep_candidates[[prep_name]], armID = armID, tot_ev = tot_ev),
        error = function(e) { errors <<- c(errors, paste0(prep_name, "/getIPD: ", conditionMessage(e))); NULL }
      )
      if (!is.null(out)) {
        tbl <- if (is.list(out) && !is.null(out$IPD)) out$IPD else out
        tbl <- data.table::as.data.table(tbl)
        data.table::setnames(tbl, tolower(names(tbl)))
        time_col <- if ("time" %in% names(tbl)) "time" else names(tbl)[1]
        status_col <- if ("status" %in% names(tbl)) "status" else if ("event" %in% names(tbl)) "event" else names(tbl)[2]
        rec <- data.table::data.table(time = as.numeric(tbl[[time_col]]), status = as.integer(tbl[[status_col]]), arm = arm)
        rec <- rec[is.finite(time) & !is.na(status)]
        if (nrow(rec) > 0L) return(rec[])
      }
    }
  }

  stop(paste(unique(errors), collapse = " | "))
}

.reconstruct_one_kmtoipd <- function(digi_dt, risk_table, arm, curve_id, censor_marks = NULL, use_marks = TRUE, drop_tol = 1e-8) {
  if (!requireNamespace("KMtoIPD", quietly = TRUE)) stop("Package KMtoIPD is not installed")
  rt <- .get_arm_risk(risk_table, arm)
  n0 <- as.integer(rt$nrisk[which.min(rt$trisk)])
  drops <- .extract_drop_bottoms(digi_dt, curve_id, drop_tol = drop_tol)
  if (nrow(drops) == 0L) {
    return(data.table::data.table(time = rep(max(rt$trisk), n0), status = 0L, arm = arm))
  }
  cens_t <- if (isTRUE(use_marks) && length(censor_marks %||% numeric(0)) > 0L) as.numeric(censor_marks) else NA
  ipd <- KMtoIPD::getIPD(n = n0, t = drops$t, S = drops$S, cens.t = cens_t)
  time <- if ("t" %in% names(ipd)) ipd$t else if ("time" %in% names(ipd)) ipd$time else ipd[[1]]
  status <- if ("event" %in% names(ipd)) ipd$event else if ("status" %in% names(ipd)) ipd$status else ipd[[2]]
  data.table::data.table(time = as.numeric(time), status = as.integer(status), arm = arm)
}

reconstruct_ipd_twoarm_final <- function(digi,
                                         risk_table,
                                         method = c("IPDfromKM", "iKM", "KMtoIPD", "QP"),
                                         curve_map = c(Control = 1, Treatment = 2),
                                         total_events = NULL,
                                         censor_marks_by_arm = NULL,
                                         kmtoipd_use_marks = TRUE,
                                         drop_tol = 1e-8) {
  method <- match.arg(method)
  digi_dt <- if (is.list(digi) && !is.null(digi$data)) digi$data else digi
  digi_dt <- .as_dt(digi_dt)
  stopifnot(all(c("time", "St", "curve") %in% names(digi_dt)))
  stopifnot(all(c("time", "nrisk_control", "nrisk_treat") %in% names(risk_table)))

  one_arm <- function(arm) {
    curve_id <- unname(curve_map[[arm]])
    marks <- if (!is.null(censor_marks_by_arm) && !is.null(censor_marks_by_arm[[arm]])) censor_marks_by_arm[[arm]] else numeric(0)
    if (method == "IPDfromKM") {
      .reconstruct_one_ipdfromkm(digi_dt, risk_table, arm, curve_id, total_events)
    } else if (method == "iKM") {
      .reconstruct_one_iKM(digi_dt, risk_table, arm, curve_id)
    } else if (method == "KMtoIPD") {
      .reconstruct_one_kmtoipd(digi_dt, risk_table, arm, curve_id, marks, use_marks = kmtoipd_use_marks, drop_tol = drop_tol)
    } else {
      .reconstruct_one_qp(digi_dt, risk_table, arm, curve_id, marks, drop_tol = drop_tol)
    }
  }

  ipd_c <- one_arm("Control")
  ipd_t <- one_arm("Treatment")
  ipd_all <- data.table::rbindlist(list(ipd_c, ipd_t), use.names = TRUE, fill = TRUE)
  ipd_all <- .fix_arm(ipd_all)
  list(method = method, ipd = ipd_all, ipd_control = ipd_c, ipd_treat = ipd_t)
}

.estimate_loghr_cox <- function(ipd) {
  ipd <- .fix_arm(ipd)
  out <- tryCatch({
    fit <- survival::coxph(survival::Surv(time, status) ~ arm, data = ipd, ties = "efron")
    cf <- stats::coef(fit)
    vc <- stats::vcov(fit)
    beta <- unname(cf["armTreatment"])
    se <- sqrt(unname(vc["armTreatment", "armTreatment"]))
    list(ok = TRUE, logHR = beta, se = se, error = NA_character_)
  }, error = function(e) list(ok = FALSE, logHR = NA_real_, se = NA_real_, error = conditionMessage(e)))
  out
}

.rmst_one <- function(ipd_arm, tau, grid_n = 1000L) {
  ipd_arm <- data.table::as.data.table(ipd_arm)
  if (nrow(ipd_arm) == 0L || !is.finite(tau) || tau <= 0) return(NA_real_)
  times <- seq(0, tau, length.out = grid_n)
  S <- .km_surv_at_times(ipd_arm, times)
  sum(diff(times) * head(S, -1L))
}

.rmst_contrast <- function(ipd, tau) {
  ipd <- .fix_arm(ipd)
  rmst_t <- .rmst_one(ipd[arm == "Treatment"], tau)
  rmst_c <- .rmst_one(ipd[arm == "Control"], tau)
  rmst_t - rmst_c
}

.km_rmse_arm <- function(true_ipd_arm, rec_ipd_arm, tau, grid_n = 200L) {
  times <- seq(0, tau, length.out = grid_n)
  S_true <- .km_surv_at_times(true_ipd_arm, times)
  S_rec <- .km_surv_at_times(rec_ipd_arm, times)
  sqrt(mean((S_rec - S_true)^2, na.rm = TRUE))
}

.digitization_rmse_arm <- function(true_ipd_arm, digi_dt, curve_id) {
  d <- .as_dt(digi_dt)[curve == curve_id]
  if (nrow(d) == 0L) return(NA_real_)
  S_true <- .km_surv_at_times(true_ipd_arm, d$time)
  sqrt(mean((d$St - S_true)^2, na.rm = TRUE))
}

.nrisk_at <- function(ipd_arm, times) {
  vapply(times, function(tt) sum(ipd_arm$time >= tt), numeric(1))
}

.surv_median <- function(ipd_arm) {
  ipd_arm <- data.table::as.data.table(ipd_arm)
  if (nrow(ipd_arm) == 0L) return(NA_real_)
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd_arm)
  med <- tryCatch(as.numeric(stats::quantile(fit, probs = 0.5)$quantile), error = function(e) NA_real_)
  med
}

evaluate_reconstruction <- function(true_ipd, rec_ipd, digi_dt, risk_table,
                                    tau_eval,
                                    curve_map = c(Control = 1, Treatment = 2)) {
  true_ipd <- .fix_arm(true_ipd)
  rec_ipd <- .fix_arm(rec_ipd)
  rt <- .as_dt(risk_table)
  tau_eval <- as.numeric(tau_eval)

  km_c <- .km_rmse_arm(true_ipd[arm == "Control"], rec_ipd[arm == "Control"], tau_eval)
  km_t <- .km_rmse_arm(true_ipd[arm == "Treatment"], rec_ipd[arm == "Treatment"], tau_eval)
  dig_c <- .digitization_rmse_arm(true_ipd[arm == "Control"], digi_dt, curve_map[["Control"]])
  dig_t <- .digitization_rmse_arm(true_ipd[arm == "Treatment"], digi_dt, curve_map[["Treatment"]])

  cox_true <- .estimate_loghr_cox(true_ipd)
  cox_rec <- .estimate_loghr_cox(rec_ipd)

  rmst_true <- .rmst_contrast(true_ipd, tau_eval)
  rmst_rec <- .rmst_contrast(rec_ipd, tau_eval)

  nrisk_ctrl_rec <- .nrisk_at(rec_ipd[arm == "Control"], rt$time)
  nrisk_trt_rec <- .nrisk_at(rec_ipd[arm == "Treatment"], rt$time)
  nrisk_ctrl_dev <- nrisk_ctrl_rec - rt$nrisk_control
  nrisk_trt_dev <- nrisk_trt_rec - rt$nrisk_treat

  true_events_c <- true_ipd[arm == "Control", sum(status == 1L)]
  true_events_t <- true_ipd[arm == "Treatment", sum(status == 1L)]
  rec_events_c <- rec_ipd[arm == "Control", sum(status == 1L)]
  rec_events_t <- rec_ipd[arm == "Treatment", sum(status == 1L)]

  med_true_c <- .surv_median(true_ipd[arm == "Control"])
  med_true_t <- .surv_median(true_ipd[arm == "Treatment"])
  med_rec_c <- .surv_median(rec_ipd[arm == "Control"])
  med_rec_t <- .surv_median(rec_ipd[arm == "Treatment"])

  data.table::data.table(
    eval_ok = isTRUE(cox_true$ok) && isTRUE(cox_rec$ok),
    digit_rmse_control = dig_c,
    digit_rmse_treatment = dig_t,
    digit_rmse_mean = mean(c(dig_c, dig_t), na.rm = TRUE),
    km_rmse_control = km_c,
    km_rmse_treatment = km_t,
    km_rmse_mean = mean(c(km_c, km_t), na.rm = TRUE),
    logHR_true = cox_true$logHR,
    logHR_rec = cox_rec$logHR,
    cox_logHR_diff = cox_rec$logHR - cox_true$logHR,
    abs_cox_logHR_diff = abs(cox_rec$logHR - cox_true$logHR),
    rmst_true = rmst_true,
    rmst_rec = rmst_rec,
    rmst_diff = rmst_rec - rmst_true,
    abs_rmst_diff = abs(rmst_rec - rmst_true),
    nrisk_control_mean_abs = mean(abs(nrisk_ctrl_dev), na.rm = TRUE),
    nrisk_treatment_mean_abs = mean(abs(nrisk_trt_dev), na.rm = TRUE),
    nrisk_mean_abs = mean(c(abs(nrisk_ctrl_dev), abs(nrisk_trt_dev)), na.rm = TRUE),
    nrisk_max_abs = max(c(abs(nrisk_ctrl_dev), abs(nrisk_trt_dev)), na.rm = TRUE),
    event_diff_control = rec_events_c - true_events_c,
    event_diff_treatment = rec_events_t - true_events_t,
    event_diff_total = (rec_events_c + rec_events_t) - (true_events_c + true_events_t),
    median_diff_control = med_rec_c - med_true_c,
    median_diff_treatment = med_rec_t - med_true_t,
    rec_n = nrow(rec_ipd),
    rec_events = rec_ipd[, sum(status == 1L)],
    rec_censors = rec_ipd[, sum(status == 0L)],
    tau_eval = tau_eval
  )
}

build_phase1_design <- function(include_stress = TRUE) {
  main <- data.table::CJ(
    scenario_id = c("W_shape1_HR067", "W_shape1_HR085", "W_shape1_HR100", "W_shape0p6_HR067", "W_shape2_HR067"),
    censoring = c("random", "exp", "front", "back", "informative"),
    sorted = FALSE
  )
  main[, `:=`(
    phase = "phase1",
    design_label = "phase1_main",
    n_label = "balanced_250_250",
    n_control = 250L,
    n_treatment = 250L,
    target_censoring = 0.30,
    digitize_sd_S = 0.01,
    digitize_mode = "event_times",
    admin_time = NA_real_,
    fixed_x_end = NA_real_
  )]

  if (!isTRUE(include_stress)) return(main[])

  stress <- data.table::CJ(
    scenario_id = c("W_shape1_HR067", "W_shape1_HR100"),
    censoring = c("random", "exp", "front", "back", "informative"),
    sorted = FALSE
  )
  stress[, `:=`(
    phase = "phase1",
    design_label = "phase1_high_censor_unbalanced",
    n_label = "unbalanced_100_200",
    n_control = 100L,
    n_treatment = 200L,
    target_censoring = 0.60,
    digitize_sd_S = 0.04,
    digitize_mode = "event_times",
    admin_time = NA_real_,
    fixed_x_end = NA_real_
  )]

  data.table::rbindlist(list(main, stress), use.names = TRUE)
}

build_phase2_design <- function(include_high_censor = TRUE) {
  main <- data.table::CJ(
    scenario_id = c("PH_H1_HR067", "PH_H2_HR085", "PH_H2_HR100", "NPHdelayed_H1_D1", "NPHdelayed_H2_D2", "NPHcross_H1_C1"),
    censoring = c("random", "front", "back"),
    sorted = FALSE
  )
  main[, `:=`(
    phase = "phase2",
    design_label = "phase2_nph_main",
    n_label = "balanced_250_250",
    n_control = 250L,
    n_treatment = 250L,
    target_censoring = 0.30,
    digitize_sd_S = 0.01,
    digitize_mode = "event_times",
    admin_time = 24,
    fixed_x_end = 24
  )]

  if (!isTRUE(include_high_censor)) return(main[])

  high <- data.table::CJ(
    scenario_id = c("NPHdelayed_H1_D1", "NPHcross_H1_C1"),
    censoring = c("random", "front", "back", "informative"),
    sorted = FALSE
  )
  high[, `:=`(
    phase = "phase2",
    design_label = "phase2_nph_high_censor",
    n_label = "balanced_250_250",
    n_control = 250L,
    n_treatment = 250L,
    target_censoring = 0.60,
    digitize_sd_S = 0.04,
    digitize_mode = "event_times",
    admin_time = 24,
    fixed_x_end = 24
  )]

  data.table::rbindlist(list(main, high), use.names = TRUE)
}

build_single_trial_design <- function(phase = c("all", "phase1", "phase2"),
                                      include_phase1_stress = TRUE,
                                      include_phase2_high_censor = TRUE) {
  phase <- match.arg(phase)
  if (phase == "phase1") return(build_phase1_design(include_stress = include_phase1_stress))
  if (phase == "phase2") return(build_phase2_design(include_high_censor = include_phase2_high_censor))
  data.table::rbindlist(list(
    build_phase1_design(include_stress = include_phase1_stress),
    build_phase2_design(include_high_censor = include_phase2_high_censor)
  ), use.names = TRUE, fill = TRUE)
}

append_fwrite <- function(dt, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(dt, file, append = file.exists(file), col.names = !file.exists(file))
}

append_log <- function(msg, file) {
  if (is.null(file) || is.na(file) || !nzchar(file)) return(invisible(NULL))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  cat(msg, "
", file = file, append = TRUE)
  invisible(NULL)
}

.with_optional_timeout <- function(expr, timeout_sec = NULL) {
  if (is.null(timeout_sec) || is.na(timeout_sec) || !is.finite(timeout_sec) || timeout_sec <= 0) {
    return(force(expr))
  }
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    warning("Package R.utils is not installed; method timeout is disabled. Install with install.packages('R.utils').", call. = FALSE)
    return(force(expr))
  }
  R.utils::withTimeout(force(expr), timeout = timeout_sec, onTimeout = "error")
}

run_one_single_trial_rep <- function(design_row, rep_id, methods = c("IPDfromKM", "iKM", "KMtoIPD", "QP"),
                                     seed_base = 20260611L,
                                     scenario_lib = get_scenario_library(),
                                     kmtoipd_use_marks = TRUE,
                                     drop_tol = 1e-8,
                                     method_timeout_sec = 180,
                                     method_log_file = NULL) {
  d <- data.table::as.data.table(design_row)
  seed <- as.integer(seed_base + as.integer(d$design_id %||% 1L) * 100000L + as.integer(rep_id))

  sim <- simulate_single_trial(
    scenario_id = d$scenario_id,
    scenario_lib = scenario_lib,
    n_control = d$n_control,
    n_treatment = d$n_treatment,
    censoring = d$censoring,
    target_censoring = d$target_censoring,
    admin_time = if (is.na(d$admin_time)) NULL else d$admin_time,
    fixed_x_end = if (is.na(d$fixed_x_end)) NULL else d$fixed_x_end,
    seed = seed
  )

  digi <- digitize_km_numeric(
    ipd_true = sim$ipd,
    risk_table = sim$risk_table,
    mode = d$digitize_mode,
    sd_S = d$digitize_sd_S,
    x_end_plot = sim$axis$x_end_plot,
    seed = seed + 10L
  )

  tau_eval <- max(sim$risk_table$time, na.rm = TRUE)

  out_list <- vector("list", length(methods))
  names(out_list) <- methods

  for (m in methods) {
    t0 <- proc.time()[[3]]
    append_log(sprintf("[%s] START design_id=%s rep=%s phase=%s scenario=%s censoring=%s label=%s method=%s",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       as.character(d$design_id %||% NA_integer_), as.character(rep_id),
                       as.character(d$phase), as.character(d$scenario_id), as.character(d$censoring),
                       as.character(d$design_label), m), method_log_file)
    rec <- tryCatch({
      .with_optional_timeout({
        reconstruct_ipd_twoarm_final(
          digi = digi,
          risk_table = sim$risk_table,
          method = m,
          curve_map = c(Control = 1, Treatment = 2),
          total_events = sim$total_events,
          censor_marks_by_arm = sim$censor_marks_by_arm,
          kmtoipd_use_marks = kmtoipd_use_marks,
          drop_tol = drop_tol
        )
      }, timeout_sec = method_timeout_sec)
    }, error = function(e) {
      list(method = m, ipd = NULL, error = conditionMessage(e))
    })
    runtime <- proc.time()[[3]] - t0
    append_log(sprintf("[%s] END   design_id=%s rep=%s method=%s runtime_sec=%.2f ok=%s error=%s",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       as.character(d$design_id %||% NA_integer_), as.character(rep_id),
                       m, as.numeric(runtime), ifelse(is.null(rec$ipd), "FALSE", "TRUE"),
                       ifelse(is.null(rec$ipd), rec$error %||% "unknown", "")), method_log_file)

    base <- data.table::copy(d)
    base[, `:=`(
      rep = as.integer(rep_id),
      seed = seed,
      method = m,
      true_n = nrow(sim$ipd),
      true_events = sim$ipd[, sum(status == 1L)],
      true_censors = sim$ipd[, sum(status == 0L)],
      achieved_censor_prop = sim$censor_prop,
      censor_scale = sim$censor_scale,
      x_end_plot = sim$axis$x_end_plot,
      runtime_sec = as.numeric(runtime)
    )]

    if (is.null(rec$ipd)) {
      base[, `:=`(ok = FALSE, error = rec$error %||% "unknown reconstruction error")]
      metrics <- data.table::data.table()
      out_list[[m]] <- base
    } else {
      met <- tryCatch({
        evaluate_reconstruction(
          true_ipd = sim$ipd,
          rec_ipd = rec$ipd,
          digi_dt = digi$data,
          risk_table = sim$risk_table,
          tau_eval = tau_eval
        )
      }, error = function(e) {
        data.table::data.table(eval_ok = FALSE, error_eval = conditionMessage(e))
      })
      base[, `:=`(ok = TRUE, error = NA_character_)]
      out_list[[m]] <- cbind(base, met)
    }
  }
  data.table::rbindlist(out_list, fill = TRUE, use.names = TRUE)
}

run_single_trial_benchmark <- function(design,
                                       n_rep = 10L,
                                       methods = c("IPDfromKM", "iKM", "KMtoIPD", "QP"),
                                       out_dir = "outputs/single_trial",
                                       results_file = file.path(out_dir, "single_trial_results.csv"),
                                       seed_base = 20260611L,
                                       overwrite = FALSE,
                                       progress_every = 10L,
                                       kmtoipd_use_marks = TRUE,
                                       drop_tol = 1e-8,
                                       method_timeout_sec = 180) {
  design <- data.table::as.data.table(design)
  if (!"design_id" %in% names(design)) design[, design_id := .I]
  data.table::setcolorder(design, c("design_id", setdiff(names(design), "design_id")))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(overwrite) && file.exists(results_file)) file.remove(results_file)
  if (isTRUE(overwrite) && file.exists(file.path(out_dir, "progress.log"))) file.remove(file.path(out_dir, "progress.log"))
  if (isTRUE(overwrite) && file.exists(file.path(out_dir, "method_progress.log"))) file.remove(file.path(out_dir, "method_progress.log"))
  data.table::fwrite(design, file.path(out_dir, "design_used.csv"))

  total <- nrow(design) * as.integer(n_rep)
  counter <- 0L
  start_time <- proc.time()[[3]]

  for (i in seq_len(nrow(design))) {
    for (r in seq_len(as.integer(n_rep))) {
      counter <- counter + 1L
      res <- tryCatch({
        run_one_single_trial_rep(
          design_row = design[i],
          rep_id = r,
          methods = methods,
          seed_base = seed_base,
          kmtoipd_use_marks = kmtoipd_use_marks,
          drop_tol = drop_tol,
          method_timeout_sec = method_timeout_sec,
          method_log_file = file.path(out_dir, "method_progress.log")
        )
      }, error = function(e) {
        bad <- data.table::copy(design[i])
        bad <- bad[rep(seq_len(.N), length(methods))]
        bad[, `:=`(
          rep = as.integer(r),
          seed = as.integer(seed_base + i * 100000L + r),
          method = methods,
          ok = FALSE,
          error = paste0("trial-level error: ", conditionMessage(e))
        )]
        bad
      })
      append_fwrite(res, results_file)
      if (counter %% progress_every == 0L || counter == total) {
        elapsed <- proc.time()[[3]] - start_time
        msg <- sprintf("[%s] Completed %d/%d design-reps (%.1f%%); elapsed %.1f min",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"), counter, total,
                       100 * counter / total, elapsed / 60)
        message(msg)
        writeLines(msg, con = file.path(out_dir, "progress.log"), sep = "\n", useBytes = TRUE)
      }
    }
  }
  invisible(data.table::fread(results_file))
}

.num_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  as.numeric(mean(x, na.rm = TRUE))
}

.num_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}

.num_q <- function(x, p = 0.90) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

.num_success <- function(x) {
  as.numeric(mean(x %in% TRUE, na.rm = TRUE))
}

summarise_single_trial_results <- function(results) {
  res <- data.table::as.data.table(results)
  if (!"ok" %in% names(res)) stop("results must contain column ok")
  if (!"eval_ok" %in% names(res)) res[, eval_ok := NA]
  group_cols <- intersect(c("phase", "design_label", "scenario_id", "n_label", "censoring", "target_censoring", "digitize_sd_S", "method"), names(res))
  out <- res[, .(
    n_runs = as.integer(.N),
    success_rate = .num_success(ok),
    eval_success_rate = .num_success(eval_ok),
    mean_digit_rmse = .num_mean(digit_rmse_mean),
    median_digit_rmse = .num_median(digit_rmse_mean),
    mean_km_rmse = .num_mean(km_rmse_mean),
    median_km_rmse = .num_median(km_rmse_mean),
    p90_km_rmse = .num_q(km_rmse_mean, 0.90),
    mean_abs_logHR_diff = .num_mean(abs_cox_logHR_diff),
    median_abs_logHR_diff = .num_median(abs_cox_logHR_diff),
    p90_abs_logHR_diff = .num_q(abs_cox_logHR_diff, 0.90),
    mean_abs_rmst_diff = .num_mean(abs_rmst_diff),
    median_abs_rmst_diff = .num_median(abs_rmst_diff),
    p90_abs_rmst_diff = .num_q(abs_rmst_diff, 0.90),
    mean_nrisk_max_abs = .num_mean(nrisk_max_abs),
    median_nrisk_max_abs = .num_median(nrisk_max_abs),
    mean_abs_event_diff_total = .num_mean(abs(event_diff_total)),
    median_abs_event_diff_total = .num_median(abs(event_diff_total)),
    median_runtime_sec = .num_median(runtime_sec),
    p90_runtime_sec = .num_q(runtime_sec, 0.90)
  ), by = group_cols]
  data.table::setorderv(out, intersect(c("phase", "design_label", "scenario_id", "method"), names(out)))
  out
}

summarise_single_trial_overall <- function(results) {
  res <- data.table::as.data.table(results)
  if (!"ok" %in% names(res)) stop("results must contain column ok")
  out <- res[, .(
    n_runs = as.integer(.N),
    success_rate = .num_success(ok),
    median_km_rmse = .num_median(km_rmse_mean),
    p90_km_rmse = .num_q(km_rmse_mean, 0.90),
    median_abs_logHR_diff = .num_median(abs_cox_logHR_diff),
    p90_abs_logHR_diff = .num_q(abs_cox_logHR_diff, 0.90),
    median_abs_rmst_diff = .num_median(abs_rmst_diff),
    p90_abs_rmst_diff = .num_q(abs_rmst_diff, 0.90),
    median_nrisk_max_abs = .num_median(nrisk_max_abs),
    median_abs_event_diff_total = .num_median(abs(event_diff_total)),
    median_runtime_sec = .num_median(runtime_sec)
  ), by = .(phase, method)]
  data.table::setorderv(out, c("phase", "method"))
  out
}

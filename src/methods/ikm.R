# reconstructKM 0.4.0 input preparation and reconstruction.
RECONSTRUCTKM_IMPLEMENTATION <- "reconstructKM::KM_reconstruct"
RECONSTRUCTKM_VERSION <- as.character(utils::packageVersion("reconstructKM"))

.reconstructkm_step_survival <- function(time, survival, query_time) {
  idx <- findInterval(query_time, time)
  idx[idx < 1L] <- 1L
  survival[pmin(idx, length(survival))]
}

.prepare_reconstructkm_direct_augmented <- function(digi_dt, risk_table, arm, curve_id) {
  xy <- data.table::as.data.table(.get_curve_xy(digi_dt, curve_id))
  data.table::setnames(xy, c("x", "y"), c("time", "surv"))
  xy <- xy[is.finite(time) & is.finite(surv)]
  xy <- xy[, .(surv = min(surv, na.rm = TRUE)), by = time]
  data.table::setorder(xy, time)
  xy[, surv := pmin(pmax(surv, 0), 1)]
  xy[, surv := cummin(surv)]

  rt <- data.table::as.data.table(.get_arm_risk(risk_table, arm))
  data.table::setnames(rt, c("trisk", "nrisk"), c("time", "NAR"))
  rt <- rt[is.finite(time) & !is.na(NAR)]
  rt <- rt[, .(NAR = as.integer(min(NAR, na.rm = TRUE))), by = time]
  data.table::setorder(rt, time)
  rt[, NAR := cummin(NAR)]

  if (nrow(rt) < 2L) stop("reconstructKM requires at least two risk-table times")
  if (rt$time[1L] != 0) stop("risk table must start at time 0")
  if (nrow(xy) == 0L) stop("no digitized coordinates")

  if (xy$time[1L] > 0 || xy$surv[1L] < 1) {
    xy <- data.table::rbindlist(list(
      data.table::data.table(time = 0, surv = 1), xy
    ))
    xy <- xy[, .(surv = min(surv, na.rm = TRUE)), by = time]
    data.table::setorder(xy, time)
    xy[, surv := cummin(surv)]
  } else {
    xy[time == 0, surv := 1]
  }

  missing_times <- setdiff(rt$time, xy$time)
  if (length(missing_times)) {
    anchors <- data.table::data.table(
      time = missing_times,
      surv = .reconstructkm_step_survival(xy$time, xy$surv, missing_times)
    )
    xy <- data.table::rbindlist(list(xy, anchors))
    xy <- xy[, .(surv = min(surv, na.rm = TRUE)), by = time]
    data.table::setorder(xy, time)
    xy[, surv := cummin(surv)]
  }

  tau <- max(rt$time)
  xy <- xy[time <= tau]

  lower <- integer(nrow(rt))
  upper <- integer(nrow(rt))
  for (ii in seq_len(nrow(rt))) {
    lower[ii] <- which(xy$time >= rt$time[ii])[1L]
    if (!is.finite(lower[ii])) stop("failed to align a risk-table time")
    if (ii < nrow(rt)) {
      candidates <- which(xy$time < rt$time[ii + 1L])
      upper[ii] <- if (length(candidates)) max(candidates) else lower[ii]
    } else {
      upper[ii] <- lower[ii]
    }
  }

  if (any(!is.finite(lower)) || any(!is.finite(upper))) {
    stop("non-finite reconstructKM interval indices")
  }
  if (any(lower > upper)) stop("invalid reconstructKM interval indices")
  if (any(diff(xy$time) <= 0)) stop("reconstructKM coordinates must have unique increasing times")

  list(
    aug_NAR = data.frame(
      time = rt$time,
      NAR = as.integer(rt$NAR),
      lower = as.integer(lower),
      upper = as.integer(upper)
    ),
    aug_surv = data.frame(time = xy$time, surv = xy$surv),
    n0 = as.integer(rt$NAR[1L]),
    tau = tau
  )
}

.prepare_reconstructkm_native <- function(digi_dt, risk_table, arm, curve_id,
                                          drop_tol = 1e-10) {
  xy <- data.table::as.data.table(.get_curve_xy(digi_dt, curve_id))
  data.table::setnames(xy, c("x", "y"), c("time", "survival"))
  xy <- xy[is.finite(time) & is.finite(survival)]
  xy <- xy[, .(survival = min(survival, na.rm = TRUE)), by = time]
  data.table::setorder(xy, time)
  xy[, survival := pmin(pmax(survival, 0), 1)]
  xy[, survival := cummin(survival)]

  rt <- data.table::as.data.table(.get_arm_risk(risk_table, arm))
  data.table::setnames(rt, c("trisk", "nrisk"), c("time", "NAR"))
  rt <- rt[is.finite(time) & !is.na(NAR)]
  rt <- rt[, .(NAR = as.integer(min(NAR, na.rm = TRUE))), by = time]
  data.table::setorder(rt, time)
  rt[, NAR := cummin(NAR)]

  if (nrow(rt) < 2L) stop("reconstructKM requires at least two risk-table times")
  if (rt$time[1L] != 0) stop("risk table must start at time 0")
  if (nrow(xy) == 0L) stop("no digitized coordinates")

  if (xy$time[1L] > 0 || xy$survival[1L] < 1) {
    xy <- data.table::rbindlist(list(
      data.table::data.table(time = 0, survival = 1), xy
    ))
    xy <- xy[, .(survival = min(survival, na.rm = TRUE)), by = time]
    data.table::setorder(xy, time)
    xy[, survival := cummin(survival)]
  } else {
    xy[time == 0, survival := 1]
  }

  tau <- max(rt$time)
  xy <- xy[time < tau]
  if (!nrow(xy)) stop("no digitized coordinates before the reconstruction horizon")
  keep <- c(TRUE, diff(xy$survival) < -abs(drop_tol))
  raw_surv <- xy[keep, .(time, survival)]
  if (nrow(raw_surv) < 2L) {
    stop("no KM decrease before the final risk-table time")
  }

  formatted <- reconstructKM::format_raw_tabs(
    raw_NAR = data.frame(time = rt$time, NAR = as.integer(rt$NAR)),
    raw_surv = data.frame(time = raw_surv$time, survival = raw_surv$survival)
  )
  list(
    aug_NAR = formatted$aug_NAR,
    aug_surv = formatted$aug_surv,
    n0 = as.integer(rt$NAR[1L]),
    tau = tau
  )
}

.validate_reconstructkm_ipd <- function(ipd, n0, tau) {
  if (!is.data.frame(ipd)) stop("reconstructKM did not return a data frame")
  if (nrow(ipd) != n0) stop("reconstructed row count differs from initial sample size")
  if (any(!is.finite(ipd$time)) || any(ipd$time < 0)) stop("invalid reconstructed follow-up times")
  if (any(!ipd$status %in% c(0L, 1L))) stop("invalid reconstructed event indicators")
  if (any(ipd$time > tau + sqrt(.Machine$double.eps))) stop("reconstructed time exceeds benchmark horizon")
  if (sum(ipd$status == 1L) > n0) stop("reconstructed event count exceeds sample size")
  invisible(TRUE)
}

.reconstruct_one_reconstructkm <- function(digi_dt, risk_table, arm, curve_id,
                                           timeout_sec = getOption("km2ipd.reconstructkm_timeout", 30)) {
  input_mode <- getOption("km2ipd.reconstructkm_input_mode", "direct_augmented")
  prepared <- if (identical(input_mode, "native")) {
    .prepare_reconstructkm_native(digi_dt, risk_table, arm, curve_id)
  } else if (identical(input_mode, "direct_augmented")) {
    .prepare_reconstructkm_direct_augmented(digi_dt, risk_table, arm, curve_id)
  } else {
    stop("Unknown reconstructKM input mode: ", input_mode)
  }
  fit <- R.utils::withTimeout(
    reconstructKM::KM_reconstruct(
      aug_NAR = prepared$aug_NAR,
      aug_surv = prepared$aug_surv
    ),
    timeout = timeout_sec,
    onTimeout = "error"
  )

  rec <- data.table::data.table(
    time = as.numeric(fit$IPD_time),
    status = as.integer(fit$IPD_event),
    arm = arm
  )
  .validate_reconstructkm_ipd(rec, prepared$n0, prepared$tau)
  data.table::setorder(rec, time, status)
  rec[]
}

.reconstruct_one_iKM <- function(digi_dt, risk_table, arm, curve_id, total_events = NULL) {
  .reconstruct_one_reconstructkm(digi_dt, risk_table, arm, curve_id)
}

#' Print an srm object
#'
#' @param x An object of class `srm`.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#'
#' @examples
#' mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
#' rownames(mat) = colnames(mat) = c("A", "B", "C")
#' fit = srm(mat)
#' print(fit)
#'
#' @export
print.srm <- function(x, ...) {
	n_time <- x$n_time
	n_act  <- x$n_actors

	cat("Social Relations Model\n")
	cat(paste0(rep("-", 40), collapse = ""), "\n")

	if (x$bipartite) {
		cat("Mode:       bipartite\n")
	} else {
		cat("Mode:       unipartite\n")
	}

	if (n_time == 1) {
		cat(sprintf("Actors:     %d\n", n_act[1]))
		cat(sprintf("Grand mean: %.4f\n", x$grand_mean[1]))
	} else {
		cat(sprintf("Time points: %d\n", n_time))
		cat(sprintf("Actors:      %s\n",
		            paste(unique(n_act), collapse = ", ")))
		cat(sprintf("Grand mean:  %.4f (avg)\n", mean(x$grand_mean)))
	}

	cat("\nVariance Decomposition:\n")

	# average across time points
	var_comps <- c("actor_var", "partner_var", "unique_var")
	var_vals <- vapply(var_comps, function(vc) {
		mean(x$stats$variance[x$stats$component == vc], na.rm = TRUE)
	}, numeric(1))

	if (any(is.na(var_vals))) {
		cat("  (Not estimable -- network too small, need n >= 4)\n")
	} else {
		has_negative <- any(var_vals < 0)
		total <- sum(pmax(var_vals, 0))
		pcts <- if (total > 0) pmax(var_vals, 0) / total * 100 else rep(NA, 3)

		labels <- c("Actor", "Partner", "Unique")
		for (i in seq_along(labels)) {
			pct_str <- if (is.na(pcts[i]) || var_vals[i] < 0) "   --" else sprintf("%5.1f%%", pcts[i])
			cat(sprintf("  %-12s %8.4f  (%s)\n", labels[i], var_vals[i], pct_str))
		}

		# covariances for unipartite
		if (!x$bipartite) {
			cov_comps <- c("relationship_cov", "actor_partner_cov")
			cov_labels <- c("Relationship", "Actor-Partner")
			for (i in seq_along(cov_comps)) {
				val <- mean(x$stats$variance[x$stats$component == cov_comps[i]],
				            na.rm = TRUE)
				cat(sprintf("  %-12s %8.4f  (cov)\n", cov_labels[i], val))
			}
		}

		if (has_negative && min(n_act) < 15) {
			cat("\nNote: Negative variance estimates can occur with small\n")
			cat("networks. Consider this exploratory.\n")
		}
	}

	invisible(x)
}


#' Summarize an srm object
#'
#' Produces a detailed variance decomposition table and, for longitudinal
#' models, per-period breakdowns.
#'
#' @param object An object of class `srm`.
#' @param ... Additional arguments (ignored).
#' @return An object of class `summary.srm`, which is a list containing:
#' \describe{
#'   \item{stats}{The full stats data frame.}
#'   \item{variance_table}{A summary table with variance, percentage, and
#'     component labels.}
#'   \item{n_time}{Number of time points.}
#'   \item{bipartite}{Logical; bipartite indicator.}
#' }
#'
#' @examples
#' mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
#' rownames(mat) = colnames(mat) = c("A", "B", "C")
#' fit = srm(mat)
#' summary(fit)
#'
#' @export
summary.srm <- function(object, ...) {
	s <- object$stats

	# average across time (replace NaN from all-NA means with NA)
	comps <- unique(s$component)
	avg_var <- vapply(comps, function(vc) {
		vals <- s$variance[s$component == vc]
		vals <- vals[!is.na(vals)]
		if (length(vals) == 0) return(NA_real_)
		mean(vals)
	}, numeric(1))

	var_comps <- c("actor_var", "partner_var", "unique_var")
	var_vals <- avg_var[var_comps]
	total <- sum(pmax(var_vals, 0), na.rm = TRUE)
	pcts <- if (!is.na(total) && total > 0) {
		pmax(var_vals, 0) / total * 100
	} else {
		rep(NA_real_, 3)
	}

	vt <- data.frame(
		component = comps,
		variance  = avg_var,
		pct       = NA_real_,
		stringsAsFactors = FALSE
	)
	vt$pct[match(var_comps, vt$component)] <- pcts
	rownames(vt) <- NULL

	out <- structure(
		list(
			stats          = s,
			variance_table = vt,
			n_time         = object$n_time,
			bipartite      = object$bipartite
		),
		class = "summary.srm"
	)
	out
}


#' Print a summary.srm object
#'
#' @param x An object of class `summary.srm`.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @export
print.summary.srm <- function(x, ...) {
	cat("Social Relations Model - Variance Decomposition\n")
	cat(paste0(rep("=", 50), collapse = ""), "\n\n")

	vt <- x$variance_table

	cat(sprintf("%-22s %10s %8s\n", "Component", "Variance", "% Total"))
	cat(paste0(rep("-", 42), collapse = ""), "\n")
	pretty <- c(
		actor_var = "Actor",
		partner_var = "Partner",
		unique_var = "Unique",
		relationship_cov = "Relationship (cov)",
		actor_partner_cov = "Actor-Partner (cov)"
	)
	if (all(is.na(vt$variance))) {
		cat("  (Not estimable -- network too small, need n >= 4)\n")
	} else {
		for (i in seq_len(nrow(vt))) {
			nm <- pretty[vt$component[i]]
			if (is.na(nm)) nm <- vt$component[i]
			pct_str <- if (!is.na(vt$pct[i]) && !is.na(vt$variance[i]) && vt$variance[i] >= 0) {
				sprintf("%5.1f%%", vt$pct[i])
			} else {
				"   --"
			}
			var_str <- if (is.na(vt$variance[i])) "        NA" else sprintf("%10.4f", vt$variance[i])
			cat(sprintf("%-22s %s %8s\n", nm, var_str, pct_str))
		}
	}

	if (x$n_time > 1) {
		cat(sprintf("\n(Averaged across %d time points)\n", x$n_time))
	}

	invisible(x)
}


#' Plot an srm object
#'
#' Dispatches to actor, partner, dyadic, or variance partition plots.
#'
#' @param x An object of class `srm`.
#' @param type One of `"actor"`, `"partner"`, `"dyadic"`, or `"variance"`.
#' @param n Maximum number of actors to display (for actor/partner/dyadic).
#' @param time Optional time point(s) to restrict to.
#' @param facet Logical; facet across time for actor/partner plots.
#' @param ... Additional arguments (ignored).
#' @return A `ggplot` object or list of `ggplot` objects.
#'
#' @examples
#' data(classroom)
#' fit = srm(classroom)
#' plot(fit, type = "actor")
#' plot(fit, type = "variance")
#'
#' @export
plot.srm <- function(x, type = c("actor", "partner", "dyadic", "variance"),
                     n = 10, time = NULL, facet = FALSE, ...) {
	type <- match.arg(type)

	if (type == "variance") {
		return(.plot_variance_partition(x, time = time))
	}

	# select the right effects
	effects <- switch(type,
		actor   = x$actor_effects,
		partner = x$partner_effects,
		dyadic  = x$unique_effects
	)

	# time filter
	if (!is.null(time)) {
		invalid_t <- time[!time %in% names(effects)]
		if (length(invalid_t) > 0) {
			cli::cli_abort("Time point(s) not found: {.val {invalid_t}}.")
		}
		effects <- effects[time]
	}

	# delegate to srm_plot internals
	srm_plot(effects, type = type, n = n, facet = facet)
}


# variance partition plot

.plot_variance_partition <- function(x, time = NULL) {
	s <- x$stats
	if (!is.null(time)) {
		s <- s[s$time %in% time, , drop = FALSE]
	}

	var_comps <- c("actor_var", "partner_var", "unique_var")
	s_var <- s[s$component %in% var_comps, , drop = FALSE]
	s_var <- s_var[!is.na(s_var$variance), , drop = FALSE]

	if (nrow(s_var) == 0) {
		cli::cli_abort("Variance components are not estimable (network too small, need n >= 4).")
	}

	# average across time
	avg <- stats::aggregate(variance ~ component, data = s_var, FUN = mean)
	avg$variance <- pmax(avg$variance, 0)
	total <- sum(avg$variance)
	avg$pct <- if (total > 0) avg$variance / total * 100 else NA

	pretty <- c(actor_var = "Actor", partner_var = "Partner",
	            unique_var = "Unique")
	avg$label <- pretty[avg$component]
	avg$label <- factor(avg$label, levels = c("Actor", "Partner", "Unique"))

	ggplot2::ggplot(avg, ggplot2::aes(x = label, y = pct)) +
		ggplot2::geom_col(width = 0.6) +
		ggplot2::geom_text(
			ggplot2::aes(label = sprintf("%.1f%%", pct)),
			vjust = -0.5, size = 3.5
		) +
		ggplot2::ylim(0, max(avg$pct, na.rm = TRUE) * 1.15) +
		ggplot2::labs(x = NULL, y = "% of Total Variance") +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			panel.grid.major.x = ggplot2::element_blank(),
			panel.grid.minor = ggplot2::element_blank()
		)
}

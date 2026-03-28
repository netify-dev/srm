#' Track SRM effects over time
#'
#' Extracts actor or partner effects across time points from a fitted
#' `srm` object and returns a tidy data frame suitable for plotting
#' temporal trends.
#'
#' @param object An object of class `srm` with multiple time points.
#' @param type One of `"actor"` or `"partner"`.
#' @param actors Optional character vector of actor names to include.
#'   If `NULL`, all actors are included.
#'
#' @return A data frame with columns: `actor`, `time`, `effect`.
#'
#' @examples
#' \donttest{
#' sim = sim_srm(n_actors = 6, n_time = 5, seed = 6886)
#' fit = srm(sim$Y)
#' trends = srm_trends(fit, type = "actor")
#' head(trends)
#' }
#'
#' @export
srm_trends <- function(object, type = c("actor", "partner"), actors = NULL) {
	if (!inherits(object, "srm")) {
		cli::cli_abort("{.arg object} must be of class {.cls srm}.")
	}
	type <- match.arg(type)

	effects <- if (type == "actor") object$actor_effects else object$partner_effects

	rows <- lapply(names(effects), function(tp) {
		e <- effects[[tp]]
		data.frame(
			actor = rownames(e),
			time  = tp,
			effect = e[, 1],
			stringsAsFactors = FALSE,
			row.names = NULL
		)
	})
	df <- do.call(rbind, rows)

	if (!is.null(actors)) {
		df <- df[df$actor %in% actors, , drop = FALSE]
	}

	df
}


#' Plot SRM effects over time
#'
#' Creates a line plot of actor or partner effects across time points.
#'
#' @param object An object of class `srm` with multiple time points.
#' @param type One of `"actor"` or `"partner"`.
#' @param actors Optional character vector of actor names to highlight.
#'   If `NULL`, all actors are shown.
#' @param n If `actors` is `NULL`, the top `n` actors by mean absolute
#'   effect are shown (default 8).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \donttest{
#' sim = sim_srm(n_actors = 8, n_time = 5, seed = 6886)
#' fit = srm(sim$Y)
#' srm_trend_plot(fit, type = "actor", n = 4)
#' }
#'
#' @importFrom tools toTitleCase
#' @export
srm_trend_plot <- function(object, type = c("actor", "partner"),
	actors = NULL, n = 8) {
	type <- match.arg(type)
	df <- srm_trends(object, type = type)

	if (is.null(actors)) {
		# select top n actors by mean absolute effect
		avg <- stats::aggregate(effect ~ actor, data = df,
			FUN = function(x) mean(abs(x)))
		avg <- avg[order(avg$effect, decreasing = TRUE), ]
		actors <- avg$actor[seq_len(min(n, nrow(avg)))]
	}

	df <- df[df$actor %in% actors, , drop = FALSE]

	ggplot2::ggplot(df, ggplot2::aes(x = time, y = effect,
		group = actor, color = actor)) +
		ggplot2::geom_line(linewidth = 0.8) +
		ggplot2::geom_point(size = 2) +
		ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
		ggplot2::labs(
			x = "Time",
			y = paste0(tools::toTitleCase(type), " Effect"),
			color = "Actor"
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			panel.grid.minor = ggplot2::element_blank()
		)
}


#' Compute stability of effects across time
#'
#' Calculates correlations of actor (or partner) effects between consecutive
#' time points, measuring how stable individual positions are over time.
#'
#' @param object An object of class `srm` with multiple time points.
#' @param type One of `"actor"` or `"partner"`.
#'
#' @return A data frame with columns: `time1`, `time2`, `correlation`, `n`.
#'   The `correlation` column is `NA` when fewer than 3 actors are shared
#'   between consecutive time points (Pearson correlation requires at
#'   least 3 observations).
#'
#' @examples
#' \donttest{
#' sim = sim_srm(n_actors = 10, n_time = 4, seed = 6886)
#' fit = srm(sim$Y)
#' srm_stability(fit, type = "actor")
#' }
#'
#' @export
srm_stability <- function(object, type = c("actor", "partner")) {
	if (!inherits(object, "srm")) {
		cli::cli_abort("{.arg object} must be of class {.cls srm}.")
	}
	type <- match.arg(type)

	effects <- if (type == "actor") object$actor_effects else object$partner_effects
	tps <- names(effects)

	if (length(tps) < 2) {
		cli::cli_abort("Need at least 2 time points to compute stability.")
	}

	rows <- list()
	for (i in seq_len(length(tps) - 1)) {
		e1 <- effects[[i]]
		e2 <- effects[[i + 1]]
		shared <- intersect(rownames(e1), rownames(e2))
		if (length(shared) < 3) {
			r <- NA_real_
		} else {
			r <- stats::cor(e1[shared, 1], e2[shared, 1])
		}
		rows[[i]] <- data.frame(
			time1 = tps[i],
			time2 = tps[i + 1],
			correlation = r,
			n = length(shared),
			stringsAsFactors = FALSE
		)
	}
	do.call(rbind, rows)
}

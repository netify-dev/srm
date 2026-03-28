#' Create plots for SRM actor, dyadic, partner effects
#'
#' Visualizes SRM effects computed by [srm_effects()]. Actor and partner
#' effects produce horizontal bar charts sorted by absolute magnitude
#' (dark bars = positive, gray bars = negative). Unique (dyadic) effects
#' produce a heatmap on a diverging color scale: blue cells indicate
#' weaker-than-expected ties, white cells are close to predicted, and
#' dark cells indicate stronger-than-expected ties.
#'
#' For the unified interface (plot directly from an `srm` object), use
#' [plot.srm()] instead.
#'
#' @param mat An effect matrix or list of effect matrices from
#'   [srm_effects()]. Must not be an `srm` object.
#' @param type One of `"actor"`, `"partner"`, or `"dyadic"`. For bipartite
#'   networks, `"dyadic"` selects top rows and columns independently.
#' @param n Maximum number of actors/partners to display, selected by
#'   absolute magnitude (default 10).
#' @param time Optional character vector of time-point names to subset
#'   when `mat` is a named list.
#' @param facet If `TRUE` and `mat` is a list, facet actor/partner plots
#'   across time. Dyadic heatmaps are returned as separate plots instead.
#' @return A `ggplot` object, or a list of `ggplot` objects for
#'   unfaceted longitudinal input.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom cli cli_abort
#' @importFrom stats reorder
#' @examples
#' # actor effect bar plot
#' data(classroom)
#' actor_eff = srm_effects(classroom, type = "actor")
#' srm_plot(actor_eff, type = "actor", n = 8)
#'
#' # dyadic heatmap
#' unique_eff = srm_effects(classroom, type = "unique")
#' srm_plot(unique_eff, type = "dyadic", n = 8)
#' @export
srm_plot <- function(mat, type = c("actor", "partner", "dyadic"), n = 10, time = NULL, facet = FALSE) {

	type <- match.arg(type)

	# validate input
	if (inherits(mat, "srm")) {
		cli::cli_abort(c(
			"Pass effect matrices from {.fn srm_effects}, not an {.cls srm} object.",
			"i" = "Use {.code plot(fit)} instead."
		))
	}
	if (!is.matrix(mat) && !is.list(mat)) {
		cli::cli_abort(c(
			"Input must be a matrix or a list of matrices.",
			"i" = "Use {.fn srm_effects} first, then pass the result to {.fn srm_plot}."
		))
	}

	# wrap single matrix in list
	if (is.matrix(mat)) {
		mat <- list(mat)
	}

	# filter by time if requested
	if (!is.null(time)) {
		if (is.null(names(mat))) {
			cli::cli_abort("Time filtering requested, but list of matrices has no names.")
		}
		invalid_t <- time[!time %in% names(mat)]
		if (length(invalid_t) > 0) {
			cli::cli_abort(paste("Time point(s) not found:", paste(invalid_t, collapse = ", ")))
		}
		mat <- mat[time]
	}

	# faceted actor/partner plot for multiple matrices
	if (facet && type %in% c("actor", "partner") && length(mat) > 1) {
		combined_data <- NULL

		for (i in seq_along(mat)) {
			x <- mat[[i]]
			period_name <- names(mat)[i]
			if (is.null(period_name)) period_name <- paste("Matrix", i)

			temp_data <- data.frame(
					id = rownames(x),
					effect = x[, 1],
					period = period_name
				)

			# sort and keep top n
			temp_data <- temp_data[order(abs(temp_data$effect), decreasing = TRUE), ]
			temp_data <- temp_data[1:min(n, nrow(temp_data)), ]

			# encode sign as shade
			temp_data$color <- ifelse(temp_data$effect > 0, "#252525", "#737373")

			if (is.null(combined_data)) {
				combined_data <- temp_data
			} else {
				combined_data <- rbind(combined_data, temp_data)
			}
		}

		p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = effect, y = reorder(id, effect), fill = color)) +
			ggplot2::geom_col() +
			ggplot2::facet_wrap(~ period, scales = "free", ncol = 2) +
			ggplot2::scale_fill_identity() +
			ggplot2::theme_bw() +
			ggplot2::theme(
				legend.position = "none",
				panel.border = ggplot2::element_blank(),
				panel.grid.major.y = ggplot2::element_blank(),
				panel.grid.major.x = ggplot2::element_blank(),
				panel.grid.minor.x = ggplot2::element_blank(),
				axis.text.y = ggplot2::element_text(size = 7),
				axis.text.x = ggplot2::element_text(size = 7, angle = 45, hjust = 1),
				plot.margin = ggplot2::unit(c(0.3, 0.5, 0.3, 0.5), "cm"),
				panel.spacing = ggplot2::unit(1, "cm"),
				strip.background = ggplot2::element_rect(fill = "black", color = "black"),
				strip.text.x = ggplot2::element_text(color = "white", hjust = 0),
				strip.text.y = ggplot2::element_text(color = "white", hjust = 0)
			) +
			ggplot2::labs(x = "Effect", y = NULL)

		return(p)
	}
	else if (facet && type == "dyadic" && length(mat) > 1) {
		cli::cli_warn("Faceting multiple dyadic heatmaps is not supported. Returning separate plots instead.")
	}

	# iterate over matrices
	res <- lapply(seq_along(mat), function(i) {

		x <- mat[[i]]

		if (type %in% c("actor", "partner")) {
			ggdata <- data.frame(
				id = rownames(x),
				effect = x[, 1]
			)

			# sort by magnitude
			ggdata <- ggdata[order(abs(ggdata$effect), decreasing = TRUE), ]

			# keep top n
			n <- min(n, nrow(ggdata))
			ggdata <- ggdata[1:n, ]

			# encode sign as shade
			ggdata$color <- ifelse(ggdata$effect > 0, "#252525", "#737373")

			plot_title <- if (length(mat) > 1) names(mat)[i] else NULL
			ggplot2::ggplot(ggdata, ggplot2::aes(x = effect, y = reorder(id, effect), fill = color)) +
				ggplot2::geom_col() +
				ggplot2::ggtitle(plot_title) +
				ggplot2::scale_fill_identity() +
				ggplot2::theme_bw() +
				ggplot2::theme(
					legend.position = "none",
					panel.border = ggplot2::element_blank(),
					panel.grid.major.y = ggplot2::element_blank(),
					panel.grid.major.x = ggplot2::element_blank(),
					panel.grid.minor.x = ggplot2::element_blank(),
					axis.text.y = ggplot2::element_text(size = 7)
				) +
				ggplot2::labs(x = "Effect", y = NULL)
		}
		else if (type == "dyadic") {
			is_bip <- nrow(x) != ncol(x)

			if (is_bip) {
				# bipartite: top row and col actors independently
				row_mag <- sort(rowSums(abs(x)), decreasing = TRUE)
				col_mag <- sort(colSums(abs(x)), decreasing = TRUE)
				top_rows <- names(row_mag)[seq_len(min(n, length(row_mag)))]
				top_cols <- names(col_mag)[seq_len(min(n, length(col_mag)))]
				x_filtered <- x[top_rows, top_cols, drop = FALSE]
			} else {
				# unipartite: top actors by combined row + col magnitude
				row_sums <- rowSums(abs(x))
				col_sums <- colSums(abs(x))
				total_effects <- row_sums + col_sums
				top_n <- min(n, length(total_effects))
				top_actors <- names(sort(total_effects, decreasing = TRUE)[1:top_n])
				x_filtered <- x[top_actors, top_actors]
			}

			# reshape for heatmap
			ggdata <- reshape2::melt(x_filtered)
			colnames(ggdata) <- c("actor", "partner", "effect")

			heat_title <- if (length(mat) > 1) names(mat)[i] else NULL
			ggplot2::ggplot(ggdata, ggplot2::aes(x = actor, y = partner, fill = effect)) +
				ggplot2::geom_tile() +
				ggplot2::ggtitle(heat_title) +
				ggplot2::scale_fill_gradient2(low = "#4393c3", mid = "white", high = "#252525", midpoint = 0) +
				ggplot2::theme_bw() +
				ggplot2::theme(
					panel.border = ggplot2::element_blank(),
					axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
					panel.grid.major = ggplot2::element_blank(),
					panel.grid.minor = ggplot2::element_blank()
				) +
				ggplot2::labs(x = "Actor", y = "Partner")
		}
	})

	if (length(res) == 1) {
		return(res[[1]])
	} else {
		return(res)
	}
}

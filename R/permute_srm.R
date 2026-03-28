#' Permutation test for SRM variance components
#'
#' Tests whether actor, partner, and unique variance components are
#' significantly different from zero using a permutation (randomization)
#' approach that preserves network size but destroys relational structure.
#'
#' @param mat A square matrix, list of square matrices, `netify` object, or
#'   an object of class `srm`. Bipartite (rectangular) networks are not
#'   supported; an error is raised if the input is non-square.
#' @param n_perms Number of permutations (default 1000).
#' @param time Optional time points to subset to.
#' @param seed Optional random seed for reproducibility.
#'
#' @return An object of class `srm_permtest` containing:
#' \describe{
#'   \item{observed}{Named vector of observed variance components.}
#'   \item{perm_dist}{Matrix of permuted variance components
#'     (`n_perms` x components).}
#'   \item{p_values}{One-tailed p-values for each component.}
#'   \item{n_perms}{Number of permutations used.}
#' }
#'
#' @details
#' For each permutation, the rows and columns of the adjacency matrix are
#' independently permuted, destroying actor and partner structure while
#' preserving the marginal distribution. The SRM variance components are
#' recomputed on each permuted matrix and the observed values are compared
#' to this null distribution.
#'
#' **Interpreting unique variance p-values:** Permutation typically
#' increases dyadic noise, so the observed unique variance is usually
#' smaller than the null distribution. A p-value of 1.0 for unique
#' variance is expected and does not imply that dyadic effects are absent.
#'
#' @examples
#' \donttest{
#' data(classroom)
#' pt = permute_srm(classroom, n_perms = 200, seed = 6886)
#' print(pt)
#' }
#'
#' @importFrom stats setNames
#' @export
permute_srm <- function(mat, n_perms = 1000, time = NULL, seed = NULL) {

	if (!is.null(seed)) set.seed(seed)

	# accept srm objects directly
	if (inherits(mat, "srm")) {
		if (mat$bipartite) {
			cli::cli_abort(c(
				"Permutation testing is only supported for unipartite (square) networks.",
				"i" = "Bipartite networks require a different permutation scheme."
			))
		}
		matrices <- mat$matrices
		# apply time filter if requested
		if (!is.null(time)) {
			invalid_t <- time[!time %in% names(matrices)]
			if (length(invalid_t) > 0) {
				cli::cli_abort("Time point(s) not found: {.val {invalid_t}}.")
			}
			matrices <- matrices[time]
		}
	} else {
		matrices <- .prepare_input(mat, time)
	}

	# reject bipartite (rectangular) matrices
	is_rect <- !all(vapply(matrices, function(m) nrow(m) == ncol(m), logical(1)))
	if (is_rect) {
		cli::cli_abort(c(
			"Permutation testing is only supported for unipartite (square) networks.",
			"i" = "Bipartite networks require a different permutation scheme."
		))
	}

	# observed
	obs_fit <- lapply(matrices, .srm_decompose)
	obs_stats <- lapply(obs_fit, function(r) {
		s <- r$stats
		stats::setNames(s$variance, s$component)
	})

	# average across time
	comps <- names(obs_stats[[1]])
	obs_avg <- vapply(comps, function(vc) {
		mean(vapply(obs_stats, function(s) s[vc], numeric(1)))
	}, numeric(1))

	# permutation loop
	perm_mat <- matrix(NA, nrow = n_perms, ncol = length(comps),
		dimnames = list(NULL, comps))

	for (p in seq_len(n_perms)) {
		perm_stats <- lapply(matrices, function(x) {
			n <- nrow(x)
			# independent row/col permutation
			row_perm <- sample.int(n)
			col_perm <- sample.int(n)
			x_perm <- x[row_perm, col_perm]
			diag(x_perm) <- 0
			rownames(x_perm) <- rownames(x)
			colnames(x_perm) <- colnames(x)
			r <- .srm_decompose(x_perm)
			stats::setNames(r$stats$variance, r$stats$component)
		})
		perm_avg <- vapply(comps, function(vc) {
			mean(vapply(perm_stats, function(s) s[vc], numeric(1)))
		}, numeric(1))
		perm_mat[p, ] <- perm_avg
	}

	# one-tailed p-values (observed >= null)
	p_values <- vapply(comps, function(vc) {
		if (is.na(obs_avg[vc])) return(NA_real_)
		(sum(perm_mat[, vc] >= obs_avg[vc], na.rm = TRUE) + 1) / (n_perms + 1)
	}, numeric(1))

	out <- structure(
		list(
			observed  = obs_avg,
			perm_dist = perm_mat,
			p_values  = p_values,
			n_perms   = n_perms
		),
		class = "srm_permtest"
	)
	out
}


#' Print a permutation test result
#'
#' Displays observed variance components, the mean of the null
#' distribution, one-tailed p-values, and significance codes.
#'
#' @param x An object of class `srm_permtest`.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @export
print.srm_permtest <- function(x, ...) {
	cat("SRM Permutation Test\n")
	cat(sprintf("Permutations: %d\n", x$n_perms))
	cat(paste0(rep("-", 50), collapse = ""), "\n")

	pretty <- c(
		actor_var = "Actor Var",
		partner_var = "Partner Var",
		unique_var = "Unique Var",
		relationship_cov = "Relationship Cov",
		actor_partner_cov = "Actor-Partner Cov"
	)

	cat(sprintf("%-22s %10s %10s %6s\n", "Component", "Observed", "Mean(Null)", "p"))
	cat(paste0(rep("-", 50), collapse = ""), "\n")

	for (vc in names(x$observed)) {
		nm <- pretty[vc]
		if (is.na(nm)) nm <- vc
		null_mean <- mean(x$perm_dist[, vc], na.rm = TRUE)
		pv <- x$p_values[vc]
		sig <- ""
		if (!is.na(pv)) {
			if (pv < 0.001) sig <- "***"
			else if (pv < 0.01) sig <- "**"
			else if (pv < 0.05) sig <- "*"
			else if (pv < 0.1) sig <- "."
		}

		obs_str <- if (is.na(x$observed[vc])) "      NA" else sprintf("%10.4f", x$observed[vc])
		null_str <- if (is.na(null_mean)) "      NA" else sprintf("%10.4f", null_mean)
		p_str <- if (is.na(pv)) "    NA" else sprintf("%6.3f", pv)

		cat(sprintf("%-22s %s %s %s %s\n", nm, obs_str, null_str, p_str, sig))
	}
	cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")

	invisible(x)
}


#' Plot permutation test null distributions
#'
#' Creates a faceted histogram of the permuted null distributions with
#' vertical dashed lines showing the observed values.
#'
#' @param x An object of class `srm_permtest`.
#' @param ... Additional arguments (ignored).
#' @return A `ggplot` object.
#' @export
plot.srm_permtest <- function(x, ...) {
	pretty <- c(
		actor_var = "Actor Var",
		partner_var = "Partner Var",
		unique_var = "Unique Var",
		relationship_cov = "Relationship Cov",
		actor_partner_cov = "Actor-Partner Cov"
	)

	# reshape perm_dist to long form
	comps <- colnames(x$perm_dist)
	long <- data.frame(
		component = rep(comps, each = x$n_perms),
		value     = as.vector(x$perm_dist),
		stringsAsFactors = FALSE
	)
	long$label <- pretty[long$component]

	obs_df <- data.frame(
		component = comps,
		value = x$observed[comps],
		label = pretty[comps],
		stringsAsFactors = FALSE
	)

	ggplot2::ggplot(long, ggplot2::aes(x = value)) +
		ggplot2::geom_histogram(bins = 30) +
		ggplot2::geom_vline(
			data = obs_df,
			ggplot2::aes(xintercept = value),
			color = "#d62728", linewidth = 0.8, linetype = "dashed"
		) +
		ggplot2::facet_wrap(~ label, scales = "free") +
		ggplot2::labs(x = "Value", y = "Count",
			title = "Permutation Null Distributions") +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0),
			panel.grid.minor = ggplot2::element_blank()
		)
}

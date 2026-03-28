#' Calculate SRM actor, partner, or unique effects
#'
#' Decomposes a round-robin network into actor (sender), partner (receiver),
#' and unique (dyadic) effects using the Social Relations Model framework.
#' Requires a square (unipartite) adjacency matrix with at least 3 actors;
#' matrices with fewer than 3 actors return `NA` values. For bipartite
#' networks, use [srm()] instead.
#'
#' @details
#' When the input matrix is symmetric (undirected), actor and partner
#' effects are identical by construction. The decomposition is most
#' informative for directed (asymmetric) networks where sender and
#' receiver roles are distinct.
#'
#' @param mat A square numeric matrix (or list of square matrices for
#'   longitudinal data). Diagonal values are set to zero internally. NAs
#'   are replaced with zero. Must not be an `srm` object (effects are
#'   already stored in `fit$actor_effects`, etc.).
#' @param type One of `"actor"`, `"partner"`, or `"unique"`:
#'   - `"actor"`: sender effects, returned as an n x 1 matrix.
#'   - `"partner"`: receiver effects, returned as an n x 1 matrix.
#'   - `"unique"`: dyadic residuals, returned as an n x n matrix with
#'     zero diagonal.
#' @param time Optional character vector of time-point names to subset
#'   when `mat` is a named list.
#' @return For a single matrix, an effect matrix (with class `actor_effect`,
#'   `partner_effect`, or `unique_effect`). For a list of matrices, a list
#'   of such matrices. See [srm()] for the unified interface.
#'
#' @examples
#' # simple inline example
#' test_matrix = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), nrow = 3, ncol = 3)
#' rownames(test_matrix) = colnames(test_matrix) = c("A", "B", "C")
#' actor_effect = srm_effects(test_matrix, type = "actor")
#'
#' # cross-sectional with included dataset
#' data(classroom)
#' actor_eff = srm_effects(classroom, type = "actor")
#' head(actor_eff)
#'
#' partner_eff = srm_effects(classroom, type = "partner")
#' unique_eff = srm_effects(classroom, type = "unique")
#'
#' # longitudinal
#' data(trade_net)
#' actor_long = srm_effects(trade_net, type = "actor")
#' actor_2017 = srm_effects(trade_net, type = "actor", time = "2017")
#'
#' @export
srm_effects <- function(mat, type = c("actor", "partner", "unique"), time = NULL) {

	# validate input
	if (inherits(mat, "srm")) {
		cli::cli_abort(c(
			"Pass a matrix or list of matrices, not an {.cls srm} object.",
			"i" = "Effects are already in {.code fit$actor_effects}, {.code fit$partner_effects}, {.code fit$unique_effects}."
		))
	}
	if (!is.matrix(mat) && !is.list(mat)) {
		cli::cli_abort(c(
			"Input must be a matrix or a list of matrices.",
			"i" = "Use {.fn netify::netify} to create a network from dyadic data."
		))
	}

	type <- match.arg(type)

	netify_warn <- FALSE

	if (is.matrix(mat)) {
		if (inherits(mat, "netify")) {
			cli::cli_alert_info("Converting 'netify' matrix to standard matrix.")
			class(mat) <- "matrix"
		}
		mat <- list(mat)
	} else if (is.list(mat)) {
		mat <- lapply(mat, function(m) {
			if (inherits(m, "netify")) {
				if (!netify_warn) {
					cli::cli_alert_info("Converting one or more 'netify' matrices in the list to standard matrices.")
					netify_warn <<- TRUE
				}
				class(m) <- "matrix"
			}
			return(m)
		})
	}

	# filter to requested time points if specified
	if (!is.null(time)) {
		if (is.null(names(mat))) {
			cli::cli_abort("Time filtering requested, but list of matrices has no names.")
		}
		invalid_t <- time[!time %in% names(mat)]
		if (length(invalid_t) > 0) {
			cli::cli_abort(paste("Time point(s) not found:",
				paste(invalid_t, collapse = ", ")))
		}
		mat <- mat[time]
	}

	results <- lapply(mat, function(x) {
		if (dim(x)[1] != dim(x)[2]) {
			cli::cli_abort(c(
				"Matrix must be square (unipartite).",
				"i" = "For bipartite (rectangular) networks, use {.fn srm} instead."
			))
		}

		# replace NAs and set diagonal to zero
		x[is.na(x)] <- 0
		diag(x) <- 0

		n <- dim(x)[1]

		# assign default names if missing
		if (is.null(rownames(x))) {
			rownames(x) <- paste0("n", sprintf("%02d", seq_len(n)))
		}
		if (is.null(colnames(x))) {
			colnames(x) <- rownames(x)
		}

		# n < 3 is too small for effect estimation
		if (n < 3) {
			if (type == "actor") {
				out <- matrix(NA_real_, n, 1,
					dimnames = list(rownames(x), "actor_effect"))
				class(out) <- c("actor_effect", class(out))
				return(out)
			} else if (type == "partner") {
				out <- matrix(NA_real_, n, 1,
					dimnames = list(rownames(x), "partner_effect"))
				class(out) <- c("partner_effect", class(out))
				return(out)
			} else {
				out <- matrix(NA_real_, n, n,
					dimnames = list(rownames(x), rownames(x)))
				if (n > 0) diag(out) <- 0
				class(out) <- c("unique_effect", class(out))
				return(out)
			}
		}

		d1 <- (n - 1)
		d2 <- (n - 2)
		d3 <- (n * d2)

		# row, column, and total means
		x_r <- rowSums(x) / d1
		x_c <- colSums(x) / d1
		x_t <- sum(x) / (n * d1)

		# actor effect: sender tendency
		a_hat <- ((d1^2) / d3 * x_r) + ((d1 / d3) * x_c) - ((d1 / d2) * x_t)

		# partner effect: receiver tendency
		b_hat <- ((d1^2) / d3 * x_c) + ((d1 / d3) * x_r) - ((d1 / d2) * x_t)

		if (type == "actor") {

			a_hat_output <- matrix(a_hat, ncol = 1,
				dimnames = list(rownames(x), "actor_effect"))
			class(a_hat_output) <- c("actor_effect", class(a_hat_output))
			return(a_hat_output)

		} else if (type == "partner") {

			b_hat_output <- matrix(b_hat, ncol = 1,
				dimnames = list(rownames(x), "partner_effect"))
			class(b_hat_output) <- c("partner_effect", class(b_hat_output))
			return(b_hat_output)

		} else if (type == "unique") {

			g_hat <- x - (
				matrix(a_hat, n, n) +
				matrix(b_hat, n, n, byrow = TRUE) +
				matrix(x_t, n, n)
			)
			diag(g_hat) <- 0
			rownames(g_hat) <- rownames(x)
			colnames(g_hat) <- rownames(x)
			class(g_hat) <- c("unique_effect", class(g_hat))
			return(g_hat)
		}

	})

	if (length(results) == 1) {
		return(results[[1]])
	} else {
		return(results)
	}
}

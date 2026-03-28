#' Calculate SRM summary statistics
#'
#' Computes descriptive statistics and variance/covariance components from a
#' round-robin network using the Social Relations Model framework. Requires
#' a square (unipartite) adjacency matrix; for bipartite networks, use
#' [srm()] instead.
#'
#' @param mat A square numeric matrix (or list of square matrices for
#'   longitudinal data). Diagonal values are set to zero internally. NAs
#'   are replaced with zero. Must not be an `srm` object (statistics are
#'   already in `summary(fit)$variance_table`).
#' @param type One of:
#'   - `"rowmeans"`: row means (average out-ties per actor).
#'   - `"colmeans"`: column means (average in-ties per actor).
#'   - `"totalmeans"`: grand mean of off-diagonal entries.
#'   - `"actor_var"`: actor (sender) variance component.
#'   - `"partner_var"`: partner (receiver) variance component.
#'   - `"unique_var"`: unique dyadic variance component.
#'   - `"relationship_cov"`: generalized reciprocity covariance.
#'   - `"actor_partner_cov"`: covariance between sending and receiving.
#'
#'   Variance and covariance components require n >= 4 actors. For n <= 3,
#'   `NA` is returned.
#' @param time Optional character vector of time-point names to subset
#'   when `mat` is a named list.
#' @return A named numeric vector (for descriptive statistics) or scalar
#'   (for variance components). For a list input, a list of such values.
#'   See [srm()] for the unified interface.
#' @importFrom cli cli_abort
#' @examples
#' # simple inline example
#' test_matrix = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), nrow = 3, ncol = 3)
#' rownames(test_matrix) = colnames(test_matrix) = c("A", "B", "C")
#' actor_totalmeans = srm_stats(test_matrix, type = "totalmeans")
#'
#' data(classroom)
#' srm_stats(classroom, type = "rowmeans")
#' srm_stats(classroom, type = "actor_var")
#' srm_stats(classroom, type = "actor_partner_cov")
#'
#' # longitudinal
#' data(trade_net)
#' srm_stats(trade_net, type = "colmeans", time = c("2015", "2019"))
#'
#' @export

srm_stats <- function(mat, type = c("rowmeans",
	"colmeans",
	"totalmeans",
	"actor_var",
	"unique_var",
	"partner_var",
	"relationship_cov",
	"actor_partner_cov"), time = NULL) {

	# validate input
	if (inherits(mat, "srm")) {
		cli::cli_abort(c(
			"Pass a matrix or list of matrices, not an {.cls srm} object.",
			"i" = "Variance components are already in {.code summary(fit)$variance_table}."
		))
	}
	if (!is.matrix(mat) && !is.list(mat)) {
		cli::cli_abort(c(
			"Input must be a matrix or a list of matrices.",
			"i" = "Use {.fn netify::netify} to create a network from dyadic data."
		))
	}

	type <- match.arg(type)

	# handle netify objects
	if (inherits(mat, "netify")) {
		cli::cli_alert_info("Converting {.cls netify} object to standard matrices.")
		if (is.matrix(mat)) {
			class(mat) <- "matrix"
		} else {
			class(mat) <- "list"
		}
	}

	if (is.matrix(mat)) {
		mat <- list(mat)
	} else if (is.list(mat)) {
		mat <- lapply(mat, function(m) {
			if (inherits(m, "netify")) class(m) <- "matrix"
			m
		})
	}

	# filter by time if specified
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

		d1 <- (n - 1)

		# n < 2 means no off-diagonal entries
		if (d1 == 0) {
			if (type %in% c("rowmeans", "colmeans")) {
				out <- rep(NA_real_, n)
				names(out) <- if (type == "rowmeans") rownames(x) else colnames(x)
				return(out)
			}
			return(NA_real_)
		}

		d2 <- (n - 2)
		d3 <- (n * d2)
		d4 <- ((d1) * (d2))
		d5 <- (((d1) * (d2) / (2)) - 1)

		x_r <- rowSums(x) / d1
		x_c <- colSums(x) / d1
		x_t <- sum(x) / (n * d1)

		if (type == "rowmeans") {
			names(x_r) <- rownames(x)
			return(x_r)
		} else if (type == "colmeans") {
			names(x_c) <- colnames(x)
			return(x_c)
		} else if (type == "totalmeans") {
			return(x_t)
		}

		a_hat <- ((d1^2) / (d3) * x_r) + (((d1) / (d3)) * x_c) - ((d1) / (d2) * x_t)
		b_hat <- ((d1^2) / (d3) * x_c) + (((d1) / (d3)) * x_r) - ((d1) / (d2) * x_t)

		g_hat <- x - (
			matrix(a_hat, n, n) +
			matrix(b_hat, n, n, byrow = TRUE) +
			matrix(x_t, n, n)
		)
		diag(g_hat) <- 0

		g_hatdiffs <- (g_hat - t(g_hat))^2
		diag(g_hatdiffs) <- 0.0
		g_hatsums <- ((g_hat + t(g_hat)) / 2)^2
		diag(g_hatsums) <- 0.0

		q1 <- sum(g_hatsums)
		q2 <- sum(g_hatdiffs) / 2

		# guard against n <= 3 where d4 or d5 is zero
		if (d4 == 0 || d5 == 0) {
			return(NA_real_)
		}

		s_2g <- ((q1 / d5) + (q2 / d4)) / 2
		s_gg <- ((q1 / d5) - (q2 / d4)) / 2

		if (type == "actor_var") {

			s_2a <- ((sum(a_hat^2) / (d1)) - ((s_2g * (d1)) / (d3)) - ((s_gg) / (d3)))
			return(s_2a)

		} else if (type == "unique_var") {

			return(s_2g)

		} else if (type == "partner_var") {

			s_2b <- ((sum(b_hat^2) / (d1)) - ((s_2g * (d1)) / (d3)) - ((s_gg) / (d3)))
			return(s_2b)

		} else if (type == "relationship_cov") {

			return(s_gg)

		} else if (type == "actor_partner_cov") {

			s_ab <- ((sum(a_hat * b_hat)) / (d1)) - ((s_gg * (d1)) / (d3)) - (s_2g / (d3))
			return(s_ab)

		}

	})

	if (length(results) == 1) {
		return(results[[1]])
	} else {
		return(results)
	}
}

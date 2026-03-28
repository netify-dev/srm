# global variables used in ggplot aes calls
utils::globalVariables(c(
	"effect", "id", "color", "actor", "partner",
	"component", "variance", "pct", "label",
	"time_point", "value", "name", "time"
))

#' Fit a Social Relations Model to network data
#'
#' Main entry point for the `srm` package. Decomposes a network (or list of
#' networks) into actor, partner, and unique dyadic effects and computes
#' the full set of SRM variance and covariance components.
#'
#' @param mat A square matrix, a list of square matrices (longitudinal),
#'   a rectangular matrix (bipartite), or a `netify` object. Diagonal
#'   values are set to zero and NAs are replaced with zero internally.
#' @param time Optional character vector of time-point names to subset to
#'   (only meaningful when `mat` is a named list).
#'
#' @details
#' The SRM decomposes each observation \eqn{X_{ij}} as:
#' \deqn{X_{ij} = \mu + a_i + b_j + g_{ij}}
#' where \eqn{\mu} is the grand mean, \eqn{a_i} is the actor (sender)
#' effect, \eqn{b_j} is the partner (receiver) effect, and \eqn{g_{ij}}
#' is the unique dyadic effect.
#'
#' For unipartite (square) networks, five variance/covariance components
#' are estimated: actor variance, partner variance, unique variance,
#' relationship covariance (generalized reciprocity), and actor-partner
#' covariance. These require at least 4 actors (n >= 4). With n = 3,
#' effects are estimated but variance components are returned as `NA`.
#'
#' Because the SRM uses method-of-moments estimation, variance estimates
#' are not constrained to be non-negative. Negative estimates can occur
#' with small networks (roughly n < 15) and indicate that the component
#' cannot be estimated reliably. When computing variance partitions,
#' negative estimates are clipped to zero before calculating percentages.
#'
#' For bipartite (rectangular) networks, only actor variance, partner
#' variance, and unique variance are estimated (covariance components
#' are not defined when senders and receivers come from different
#' populations). Effects are simple deviations from the grand mean.
#' Variance components are bias-corrected using two-way ANOVA
#' degrees of freedom. There is no minimum size restriction beyond
#' n >= 2 in each mode.
#'
#' **Symmetric networks:** When the input matrix is symmetric (undirected),
#' actor and partner effects are identical by construction. The variance
#' decomposition is constrained accordingly. The SRM is most informative
#' for directed (asymmetric) networks.
#'
#' @return An object of class `srm`, which is a list containing:
#' \describe{
#'   \item{actor_effects}{Actor (sender) effects for each time point.}
#'   \item{partner_effects}{Partner (receiver) effects for each time point.}
#'   \item{unique_effects}{Unique dyadic effects for each time point.}
#'   \item{stats}{A data frame of variance / covariance components for each
#'     time point.}
#'   \item{grand_mean}{The overall network mean for each time point.}
#'   \item{matrices}{The (possibly time-filtered) input matrices used.}
#'   \item{n_time}{Number of time points.}
#'   \item{n_actors}{Number of actors (per time point).}
#'   \item{bipartite}{Logical; whether the model was fit as bipartite.}
#' }
#'
#' @seealso [srm_effects()] for effect extraction, [srm_stats()] for
#'   individual statistics, [permute_srm()] for inference,
#'   [srm_trends()] and [srm_stability()] for longitudinal analysis,
#'   [sim_srm()] for simulation.
#'
#' @importFrom netify netify
#' @examples
#' mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
#' rownames(mat) = colnames(mat) = c("A", "B", "C")
#' fit = srm(mat)
#' fit
#'
#' data(classroom)
#' fit = srm(classroom)
#' fit
#' summary(fit)
#'
#' @export
srm <- function(mat, time = NULL) {

	# handle input
	mat <- .prepare_input(mat, time)

	# detect bipartite
	bipartite <- !all(vapply(mat, function(m) nrow(m) == ncol(m), logical(1)))

	# compute per-period results
	res <- lapply(mat, function(x) {
		if (bipartite) {
			.srm_bipartite_decompose(x)
		} else {
			.srm_decompose(x)
		}
	})

	# assemble output
	actor_effects  <- lapply(res, `[[`, "actor")
	partner_effects <- lapply(res, `[[`, "partner")
	unique_effects <- lapply(res, `[[`, "unique")
	grand_mean     <- vapply(res, `[[`, numeric(1), "grand_mean")

	# build stats data frame
	stats_list <- lapply(seq_along(res), function(i) {
		s <- res[[i]]$stats
		s$time <- names(mat)[i]
		s
	})
	stats_df <- do.call(rbind, stats_list)

	n_actors <- vapply(mat, nrow, integer(1))

	out <- structure(
		list(
			actor_effects   = actor_effects,
			partner_effects = partner_effects,
			unique_effects  = unique_effects,
			stats           = stats_df,
			grand_mean      = grand_mean,
			matrices        = mat,
			n_time          = length(mat),
			n_actors        = n_actors,
			bipartite       = bipartite
		),
		class = "srm"
	)
	out
}


# standard (square) decomposition

.srm_decompose <- function(x) {
	# assign default names if missing
	n <- nrow(x)
	if (is.null(rownames(x))) {
		rownames(x) <- paste0("n", sprintf("%02d", seq_len(n)))
	}
	if (is.null(colnames(x))) {
		colnames(x) <- rownames(x)
	}

	# n < 3 is too small for SRM effect estimation
	if (n < 3) {
		a_mat <- matrix(NA_real_, nrow = n, ncol = 1,
			dimnames = list(rownames(x), "actor_effect"))
		b_mat <- matrix(NA_real_, nrow = n, ncol = 1,
			dimnames = list(rownames(x), "partner_effect"))
		g_hat <- matrix(NA_real_, n, n,
			dimnames = list(rownames(x), rownames(x)))
		if (n > 0) diag(g_hat) <- 0
		stats_df <- data.frame(
			component = c("actor_var", "partner_var", "unique_var",
				"relationship_cov", "actor_partner_cov"),
			variance = rep(NA_real_, 5),
			stringsAsFactors = FALSE
		)
		x_t <- if (n > 1) { diag(x) <- 0; sum(x) / (n * (n - 1)) } else NA_real_
		return(list(
			actor = a_mat, partner = b_mat, unique = g_hat,
			grand_mean = x_t, stats = stats_df
		))
	}

	diag(x) <- 0
	d1 <- n - 1
	d2 <- n - 2
	d3 <- n * d2

	x_r <- rowSums(x) / d1
	x_c <- colSums(x) / d1
	x_t <- sum(x) / (n * d1)

	a_hat <- ((d1^2) / d3) * x_r + (d1 / d3) * x_c - (d1 / d2) * x_t
	b_hat <- ((d1^2) / d3) * x_c + (d1 / d3) * x_r - (d1 / d2) * x_t

	g_hat <- x - (
		matrix(a_hat, n, n) +
		matrix(b_hat, n, n, byrow = TRUE) +
		matrix(x_t, n, n)
	)
	diag(g_hat) <- 0

	# variance components
	d4 <- d1 * d2
	d5 <- d1 * d2 / 2 - 1

	g_hatdiffs <- (g_hat - t(g_hat))^2
	diag(g_hatdiffs) <- 0
	g_hatsums <- ((g_hat + t(g_hat)) / 2)^2
	diag(g_hatsums) <- 0

	q1 <- sum(g_hatsums)
	q2 <- sum(g_hatdiffs) / 2

	# guard against division by zero when n <= 3
	if (d4 == 0 || d5 == 0) {
		s_2g <- NA_real_
		s_gg <- NA_real_
		s_2a <- NA_real_
		s_2b <- NA_real_
		s_ab <- NA_real_
	} else {
		s_2g <- (q1 / d5 + q2 / d4) / 2
		s_gg <- (q1 / d5 - q2 / d4) / 2

		s_2a <- sum(a_hat^2) / d1 - s_2g * d1 / d3 - s_gg / d3
		s_2b <- sum(b_hat^2) / d1 - s_2g * d1 / d3 - s_gg / d3
		s_ab <- sum(a_hat * b_hat) / d1 - s_gg * d1 / d3 - s_2g / d3
	}

	# name effects
	a_mat <- matrix(a_hat, ncol = 1, dimnames = list(rownames(x), "actor_effect"))
	b_mat <- matrix(b_hat, ncol = 1, dimnames = list(rownames(x), "partner_effect"))
	rownames(g_hat) <- rownames(x)
	colnames(g_hat) <- rownames(x)

	stats_df <- data.frame(
		component = c("actor_var", "partner_var", "unique_var",
		              "relationship_cov", "actor_partner_cov"),
		variance  = c(s_2a, s_2b, s_2g, s_gg, s_ab),
		stringsAsFactors = FALSE
	)

	list(
		actor      = a_mat,
		partner    = b_mat,
		unique     = g_hat,
		grand_mean = x_t,
		stats      = stats_df
	)
}


# bipartite decomposition

.srm_bipartite_decompose <- function(x) {
	nr <- nrow(x)
	nc <- ncol(x)
	if (is.null(rownames(x))) {
		rownames(x) <- paste0("r", sprintf("%02d", seq_len(nr)))
	}
	if (is.null(colnames(x))) {
		colnames(x) <- paste0("c", sprintf("%02d", seq_len(nc)))
	}

	x_r <- rowMeans(x, na.rm = TRUE)
	x_c <- colMeans(x, na.rm = TRUE)
	x_t <- mean(x, na.rm = TRUE)

	a_hat <- x_r - x_t
	b_hat <- x_c - x_t

	g_hat <- x - (
		matrix(a_hat, nr, nc) +
		matrix(b_hat, nr, nc, byrow = TRUE) +
		matrix(x_t, nr, nc)
	)

	# bias-corrected variance components (two-way anova)
	s_2g <- sum(g_hat^2) / ((nr - 1) * (nc - 1))
	s_2a <- stats::var(a_hat) - s_2g / nc
	s_2b <- stats::var(b_hat) - s_2g / nr

	a_mat <- matrix(a_hat, ncol = 1,
	                dimnames = list(rownames(x), "actor_effect"))
	b_mat <- matrix(b_hat, ncol = 1,
	                dimnames = list(colnames(x), "partner_effect"))
	rownames(g_hat) <- rownames(x)
	colnames(g_hat) <- colnames(x)

	stats_df <- data.frame(
		component = c("actor_var", "partner_var", "unique_var"),
		variance  = c(s_2a, s_2b, s_2g),
		stringsAsFactors = FALSE
	)

	list(
		actor      = a_mat,
		partner    = b_mat,
		unique     = g_hat,
		grand_mean = x_t,
		stats      = stats_df
	)
}


# input preparation

.prepare_input <- function(mat, time = NULL) {

	# handle netify objects
	if (inherits(mat, "netify")) {
		cli::cli_alert_info("Converting {.cls netify} object to standard matrices.")
		mat <- .netify_to_matrices(mat)
	}

	if (is.matrix(mat)) {
		mat <- list(mat)
	}

	if (!is.list(mat) || length(mat) == 0) {
		cli::cli_abort(c(
			"{.arg mat} must be a matrix, list of matrices, or a {.cls netify} object.",
			"i" = "Use {.fn netify::netify} to create a network from dyadic data."
		))
	}

	# convert any remaining netify matrices inside the list
	na_warned <- FALSE
	mat <- lapply(mat, function(m) {
		if (inherits(m, "netify")) {
			class(m) <- setdiff(class(m), "netify")
		}
		if (!is.matrix(m)) {
			cli::cli_abort("Each element of {.arg mat} must be a matrix.")
		}
		# replace NAs with zero (e.g. from netify diag_to_NA)
		if (anyNA(m)) {
			if (!na_warned) {
				cli::cli_alert_info("NAs found in input; replacing with zero.")
				na_warned <<- TRUE
			}
			m[is.na(m)] <- 0
		}
		m
	})

	# time filtering
	if (!is.null(time)) {
		if (is.null(names(mat))) {
			cli::cli_abort("Time filtering requested but the list has no names.")
		}
		invalid_t <- time[!time %in% names(mat)]
		if (length(invalid_t) > 0) {
			cli::cli_abort("Time point(s) not found: {.val {invalid_t}}.")
		}
		mat <- mat[time]
	}

	# give default names if unnamed
	if (is.null(names(mat))) {
		names(mat) <- paste0("t", seq_along(mat))
	}

	mat
}


# netify to matrix list

.netify_to_matrices <- function(net) {
	# netify cross-sectional objects are matrices with class "netify"
	# netify longitudinal objects are lists of such matrices
	if (is.matrix(net)) {
		m <- net
		class(m) <- "matrix"
		return(list(m))
	}
	if (is.list(net)) {
		return(lapply(net, function(m) {
			if (inherits(m, "netify")) class(m) <- "matrix"
			m
		}))
	}
	cli::cli_abort("Cannot convert this {.cls netify} object to matrices.")
}

#' Simulate network data from a Social Relations Model
#'
#' Generates synthetic network matrices from known actor, partner, and unique
#' effects, useful for teaching, validation, and power analysis.
#'
#' @param n_actors Number of actors (nodes) in the network. For bipartite
#'   networks, a length-2 vector `c(n_row, n_col)` specifying the number
#'   of senders and receivers.
#' @param n_time Number of time periods (default 1 for cross-sectional).
#' @param actor_var Variance of actor (sender) effects.
#' @param partner_var Variance of partner (receiver) effects.
#' @param unique_var Variance of unique dyadic effects.
#' @param actor_partner_cov Covariance between actor and partner effects.
#'   Must satisfy `actor_partner_cov^2 <= actor_var * partner_var`.
#' @param relationship_cov Covariance of unique effects within dyads
#'   (generalized reciprocity). Must satisfy
#'   `abs(relationship_cov) <= unique_var`.
#' @param grand_mean Overall network mean.
#' @param bipartite If `TRUE`, generate a bipartite (two-mode) network.
#'   When bipartite, `n_actors` is a length-2 vector `c(n_row, n_col)`.
#'   The `actor_partner_cov` and `relationship_cov` parameters are
#'   ignored for bipartite networks (these quantities are not defined
#'   when senders and receivers are from different populations).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y}{A matrix or named list of matrices (if `n_time > 1`).}
#'   \item{truth}{A list containing the true `actor_effects`,
#'     `partner_effects`, `unique_effects`, and `grand_mean`.}
#'   \item{params}{The parameter values used for simulation.}
#' }
#'
#' @examples
#' # cross-sectional
#' sim = sim_srm(n_actors = 10, actor_var = 1, partner_var = 0.5,
#'                unique_var = 2, seed = 6886)
#' fit = srm(sim$Y)
#' summary(fit)
#'
#' \donttest{
#' # longitudinal
#' sim_long = sim_srm(n_actors = 8, n_time = 5, seed = 6886)
#' fit_long = srm(sim_long$Y)
#' fit_long
#'
#' # bipartite
#' sim_bip = sim_srm(n_actors = c(6, 8), bipartite = TRUE, seed = 6886)
#' fit_bip = srm(sim_bip$Y)
#' summary(fit_bip)
#' }
#'
#' @importFrom stats rnorm
#' @export
sim_srm <- function(n_actors = 10,
	n_time = 1,
	actor_var = 1,
	partner_var = 1,
	unique_var = 1,
	actor_partner_cov = 0,
	relationship_cov = 0,
	grand_mean = 0,
	bipartite = FALSE,
	seed = NULL) {

	if (!is.null(seed)) set.seed(seed)

	if (bipartite) {
		return(.sim_bipartite(n_actors, n_time, actor_var, partner_var,
			unique_var, grand_mean))
	}

	n <- n_actors
	if (length(n) != 1 || n < 3) {
		cli::cli_abort("{.arg n_actors} must be a single integer >= 3.")
	}

	# correlated actor/partner draws
	sigma_ab <- matrix(c(actor_var, actor_partner_cov,
		actor_partner_cov, partner_var), 2, 2)
	if (all(sigma_ab == 0)) {
		L <- matrix(0, 2, 2)
	} else {
		L <- tryCatch(chol(sigma_ab), error = function(e) {
			cli::cli_abort(c(
				"Actor-partner covariance matrix is not positive semi-definite.",
				"i" = "Ensure that {.arg actor_partner_cov}^2 <= {.arg actor_var} * {.arg partner_var}."
			))
		})
	}

	# correlated unique effects covariance
	if (abs(relationship_cov) > unique_var) {
		cli::cli_abort(c(
			"Unique effects covariance matrix is not positive semi-definite.",
			"i" = "Ensure that abs({.arg relationship_cov}) <= {.arg unique_var}."
		))
	}
	sigma_g <- matrix(c(unique_var, relationship_cov,
		relationship_cov, unique_var), 2, 2)
	if (all(sigma_g == 0)) {
		L_g <- matrix(0, 2, 2)
	} else {
		# eigendecomposition handles both PD and PSD cases
		eg <- eigen(sigma_g, symmetric = TRUE)
		eg$values <- pmax(eg$values, 0)
		L_g <- diag(sqrt(eg$values)) %*% t(eg$vectors)
	}

	sim_one <- function(period_name) {
		z <- matrix(stats::rnorm(n * 2), n, 2)
		ab <- z %*% L
		a <- ab[, 1]
		b <- ab[, 2]

		# correlated unique effects
		g <- matrix(0, n, n)
		for (i in 1:(n - 1)) {
			for (j in (i + 1):n) {
				pair <- stats::rnorm(2) %*% L_g
				g[i, j] <- pair[1]
				g[j, i] <- pair[2]
			}
		}

		Y <- grand_mean +
			matrix(a, n, n) +
			matrix(b, n, n, byrow = TRUE) +
			g
		diag(Y) <- 0

		actor_names <- paste0("n", sprintf("%02d", seq_len(n)))
		rownames(Y) <- colnames(Y) <- actor_names

		list(
			Y = Y,
			actor_effects = stats::setNames(a, actor_names),
			partner_effects = stats::setNames(b, actor_names),
			unique_effects = g,
			grand_mean = grand_mean
		)
	}

	periods <- if (n_time > 1) {
		paste0("t", seq_len(n_time))
	} else {
		"t1"
	}

	sims <- lapply(periods, sim_one)

	if (n_time == 1) {
		Y <- sims[[1]]$Y
		truth <- list(
			actor_effects   = sims[[1]]$actor_effects,
			partner_effects = sims[[1]]$partner_effects,
			unique_effects  = sims[[1]]$unique_effects,
			grand_mean      = grand_mean
		)
	} else {
		Y <- stats::setNames(lapply(sims, `[[`, "Y"), periods)
		truth <- list(
			actor_effects   = lapply(sims, `[[`, "actor_effects"),
			partner_effects = lapply(sims, `[[`, "partner_effects"),
			unique_effects  = lapply(sims, `[[`, "unique_effects"),
			grand_mean      = grand_mean
		)
	}

	list(
		Y = Y,
		truth = truth,
		params = list(
			n_actors = n_actors,
			n_time = n_time,
			actor_var = actor_var,
			partner_var = partner_var,
			unique_var = unique_var,
			actor_partner_cov = actor_partner_cov,
			relationship_cov = relationship_cov,
			grand_mean = grand_mean,
			bipartite = FALSE
		)
	)
}


# bipartite simulation helper

.sim_bipartite <- function(n_actors, n_time, actor_var, partner_var,
	unique_var, grand_mean) {
	if (length(n_actors) != 2) {
		cli::cli_abort(
			"For bipartite simulation, {.arg n_actors} must be a length-2 vector c(n_row, n_col)."
		)
	}
	nr <- n_actors[1]
	nc <- n_actors[2]

	sim_one <- function(period_name) {
		a <- stats::rnorm(nr, sd = sqrt(actor_var))
		b <- stats::rnorm(nc, sd = sqrt(partner_var))
		g <- matrix(stats::rnorm(nr * nc, sd = sqrt(unique_var)), nr, nc)

		Y <- grand_mean +
			matrix(a, nr, nc) +
			matrix(b, nr, nc, byrow = TRUE) +
			g

		rn <- paste0("r", sprintf("%02d", seq_len(nr)))
		cn <- paste0("c", sprintf("%02d", seq_len(nc)))
		rownames(Y) <- rn
		colnames(Y) <- cn

		list(
			Y = Y,
			actor_effects = stats::setNames(a, rn),
			partner_effects = stats::setNames(b, cn),
			unique_effects = g,
			grand_mean = grand_mean
		)
	}

	periods <- if (n_time > 1) paste0("t", seq_len(n_time)) else "t1"
	sims <- lapply(periods, sim_one)

	if (n_time == 1) {
		Y <- sims[[1]]$Y
		truth <- list(
			actor_effects   = sims[[1]]$actor_effects,
			partner_effects = sims[[1]]$partner_effects,
			unique_effects  = sims[[1]]$unique_effects,
			grand_mean      = grand_mean
		)
	} else {
		Y <- stats::setNames(lapply(sims, `[[`, "Y"), periods)
		truth <- list(
			actor_effects   = lapply(sims, `[[`, "actor_effects"),
			partner_effects = lapply(sims, `[[`, "partner_effects"),
			unique_effects  = lapply(sims, `[[`, "unique_effects"),
			grand_mean      = grand_mean
		)
	}

	list(
		Y = Y,
		truth = truth,
		params = list(
			n_actors = n_actors,
			n_time = n_time,
			actor_var = actor_var,
			partner_var = partner_var,
			unique_var = unique_var,
			grand_mean = grand_mean,
			bipartite = TRUE
		)
	)
}

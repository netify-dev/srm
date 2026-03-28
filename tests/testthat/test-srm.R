test_that("srm() works on a single matrix", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	fit = srm(mat)

	expect_s3_class(fit, "srm")
	expect_equal(fit$n_time, 1)
	expect_equal(fit$n_actors, c(t1 = 3L))
	expect_false(fit$bipartite)

	# check all components exist

	expect_true(all(c("actor_effects", "partner_effects", "unique_effects",
	                   "stats", "grand_mean", "matrices") %in% names(fit)))
})

test_that("srm() works on a list of matrices", {
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = matrix(c(0, 3, 1, 1, 0, 2, 2, 1, 0), 3, 3)
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	fit = srm(list("t1" = m1, "t2" = m2))

	expect_s3_class(fit, "srm")
	expect_equal(fit$n_time, 2)
	expect_equal(nrow(fit$stats), 10) # 5 components x 2 time points
})

test_that("srm() time filtering works", {
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = matrix(c(0, 3, 1, 1, 0, 2, 2, 1, 0), 3, 3)
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	fit = srm(list("2020" = m1, "2021" = m2), time = "2021")

	expect_equal(fit$n_time, 1)
	expect_equal(names(fit$matrices), "2021")
})

test_that("srm() rejects non-matrix input", {
	expect_error(srm("not a matrix"), class = "rlang_error")
	expect_error(srm(42), class = "rlang_error")
})

test_that("srm() grand mean is correct", {
	mat = matrix(c(0, 4, 6, 2, 0, 8, 10, 12, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	fit = srm(mat)
	expected_mean = sum(mat) / (3 * 2)  # n * (n-1)
	expect_equal(unname(fit$grand_mean[1]), expected_mean, tolerance = 1e-10)
})

test_that("srm() decomposition sums back to original", {
	mat = matrix(c(0, 3, 5, 7, 0, 2, 4, 6, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	fit = srm(mat)

	a = fit$actor_effects[[1]][, 1]
	b = fit$partner_effects[[1]][, 1]
	g = fit$unique_effects[[1]]
	mu = fit$grand_mean[1]

	n = 3
	reconstructed = mu +
		matrix(a, n, n) +
		matrix(b, n, n, byrow = TRUE) +
		g

	# off-diagonals should match
	for (i in 1:n) {
		for (j in 1:n) {
			if (i != j) {
				expect_equal(reconstructed[i, j], mat[i, j], tolerance = 1e-10)
			}
		}
	}
})

test_that("srm() effects sum to approximately zero", {
	sim = sim_srm(n_actors = 12, seed = 6886)
	fit = srm(sim$Y)

	a = fit$actor_effects[[1]][, 1]
	b = fit$partner_effects[[1]][, 1]

	expect_true(abs(sum(a)) < 1e-10)
	expect_true(abs(sum(b)) < 1e-10)
})

test_that("symmetric matrix gives identical actor and partner effects", {
	set.seed(6886)
	n = 8
	mat = matrix(rnorm(n * n), n, n)
	mat = (mat + t(mat)) / 2
	diag(mat) = 0
	rownames(mat) = colnames(mat) = paste0("n", 1:n)

	fit = srm(mat)
	a = fit$actor_effects[[1]][, 1]
	b = fit$partner_effects[[1]][, 1]

	expect_equal(a, b, tolerance = 1e-10)
})

test_that("hand-calculated 4x4 variance components match", {
	# verify formulas with a small hand-traceable example
	mat = matrix(c(0, 3, 5, 2,
	               7, 0, 4, 6,
	               1, 8, 0, 3,
	               4, 2, 7, 0), 4, 4, byrow = TRUE)
	rownames(mat) = colnames(mat) = c("A", "B", "C", "D")

	fit = srm(mat)

	# manually compute expected grand mean
	n = 4
	d1 = n - 1
	x_t = sum(mat) / (n * d1)
	expect_equal(unname(fit$grand_mean[1]), x_t)

	# manually compute actor effects for first actor
	x_r = rowSums(mat) / d1
	x_c = colSums(mat) / d1
	d2 = n - 2
	d3 = n * d2
	a_hat = ((d1^2) / d3) * x_r + (d1 / d3) * x_c - (d1 / d2) * x_t
	expect_equal(fit$actor_effects[[1]][, 1], a_hat, tolerance = 1e-10)

	# variance components should be finite
	expect_true(all(is.finite(fit$stats$variance)))
})

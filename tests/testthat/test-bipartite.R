test_that("srm() handles bipartite matrices", {
	set.seed(6886)
	mat = matrix(rnorm(20), 4, 5)
	rownames(mat) = paste0("S", 1:4)
	colnames(mat) = paste0("R", 1:5)

	fit = srm(mat)

	expect_s3_class(fit, "srm")
	expect_true(fit$bipartite)
	expect_equal(nrow(fit$actor_effects[[1]]), 4)
	expect_equal(nrow(fit$partner_effects[[1]]), 5)
	expect_equal(dim(fit$unique_effects[[1]]), c(4, 5))
})

test_that("bipartite variance components are computed", {
	set.seed(6886)
	mat = matrix(rnorm(20), 4, 5)
	rownames(mat) = paste0("S", 1:4)
	colnames(mat) = paste0("R", 1:5)

	fit = srm(mat)
	s = summary(fit)

	# bipartite has 3 components only
	expect_equal(nrow(s$variance_table), 3)
})

test_that("bipartite print works", {
	set.seed(6886)
	mat = matrix(rnorm(20), 4, 5)
	rownames(mat) = paste0("S", 1:4)
	colnames(mat) = paste0("R", 1:5)

	fit = srm(mat)
	out = capture.output(print(fit))
	expect_true(any(grepl("bipartite", out)))
})

test_that("bipartite decomposition sums back to original", {
	skip_on_cran()
	set.seed(6886)
	mat = matrix(rnorm(20), 4, 5)
	rownames(mat) = paste0("S", 1:4)
	colnames(mat) = paste0("R", 1:5)

	fit = srm(mat)
	a = fit$actor_effects[[1]][, 1]
	b = fit$partner_effects[[1]][, 1]
	g = fit$unique_effects[[1]]
	mu = fit$grand_mean[1]

	nr = nrow(mat)
	nc = ncol(mat)
	reconstructed = mu +
		matrix(a, nr, nc) +
		matrix(b, nr, nc, byrow = TRUE) +
		g

	for (i in 1:nr) {
		for (j in 1:nc) {
			expect_equal(reconstructed[i, j], mat[i, j], tolerance = 1e-10)
		}
	}
})

test_that("sim_srm bipartite round-trips through srm()", {
	sim = sim_srm(n_actors = c(6, 8), bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)

	expect_true(fit$bipartite)
	expect_equal(nrow(fit$actor_effects[[1]]), 6)
	expect_equal(nrow(fit$partner_effects[[1]]), 8)
})

test_that("bipartite longitudinal works", {
	skip_on_cran()
	sim = sim_srm(n_actors = c(5, 7), n_time = 3, bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)

	expect_equal(fit$n_time, 3)
	expect_true(fit$bipartite)
})

test_that("bipartite variance uses correct degrees of freedom", {
	# verify the bias-corrected estimator against known formulas
	set.seed(6886)
	nr = 6
	nc = 8
	mat = matrix(rnorm(nr * nc), nr, nc)
	rownames(mat) = paste0("S", 1:nr)
	colnames(mat) = paste0("R", 1:nc)

	fit = srm(mat)
	s = fit$stats

	# manually compute expected unique variance using correct df
	x_r = rowMeans(mat)
	x_c = colMeans(mat)
	x_t = mean(mat)
	a_hat = x_r - x_t
	b_hat = x_c - x_t
	g_hat = mat - matrix(a_hat, nr, nc) - matrix(b_hat, nr, nc, byrow = TRUE) - x_t
	expected_s2g = sum(g_hat^2) / ((nr - 1) * (nc - 1))

	actual_s2g = s$variance[s$component == "unique_var"]
	expect_equal(actual_s2g, expected_s2g, tolerance = 1e-10)

	# actor variance should subtract noise bias
	expected_s2a = var(a_hat) - expected_s2g / nc
	actual_s2a = s$variance[s$component == "actor_var"]
	expect_equal(actual_s2a, expected_s2a, tolerance = 1e-10)

	# partner variance should subtract noise bias
	expected_s2b = var(b_hat) - expected_s2g / nr
	actual_s2b = s$variance[s$component == "partner_var"]
	expect_equal(actual_s2b, expected_s2b, tolerance = 1e-10)
})

test_that("bipartite variance recovery is approximately unbiased", {
	skip_on_cran()
	# average over many replications to check for systematic bias
	set.seed(6886)
	n_reps = 200
	est_a = est_b = est_g = numeric(n_reps)
	for (i in 1:n_reps) {
		sim_i = sim_srm(n_actors = c(15, 20), bipartite = TRUE,
		                actor_var = 2.0, partner_var = 1.0, unique_var = 1.5)
		fit_i = srm(sim_i$Y)
		vt = fit_i$stats
		est_a[i] = vt$variance[vt$component == "actor_var"]
		est_b[i] = vt$variance[vt$component == "partner_var"]
		est_g[i] = vt$variance[vt$component == "unique_var"]
	}

	# mean estimates should be within 15% of true values
	expect_true(abs(mean(est_a) - 2.0) / 2.0 < 0.15)
	expect_true(abs(mean(est_b) - 1.0) / 1.0 < 0.15)
	expect_true(abs(mean(est_g) - 1.5) / 1.5 < 0.15)
})

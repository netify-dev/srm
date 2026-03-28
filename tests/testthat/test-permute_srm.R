test_that("permute_srm returns correct structure", {
	skip_on_cran()
	sim = sim_srm(n_actors = 8, seed = 6886)
	pt = permute_srm(sim$Y, n_perms = 50, seed = 6886)

	expect_s3_class(pt, "srm_permtest")
	expect_equal(pt$n_perms, 50)
	expect_equal(length(pt$observed), 5)
	expect_equal(nrow(pt$perm_dist), 50)
	expect_equal(ncol(pt$perm_dist), 5)
	expect_equal(length(pt$p_values), 5)
})

test_that("permute_srm p-values are in [0, 1]", {
	skip_on_cran()
	sim = sim_srm(n_actors = 8, seed = 6886)
	pt = permute_srm(sim$Y, n_perms = 50, seed = 6886)

	expect_true(all(pt$p_values >= 0 & pt$p_values <= 1))
})

test_that("permute_srm accepts srm objects", {
	skip_on_cran()
	sim = sim_srm(n_actors = 8, seed = 6886)
	fit = srm(sim$Y)
	pt = permute_srm(fit, n_perms = 30, seed = 6886)

	expect_s3_class(pt, "srm_permtest")
})

test_that("permute_srm print works", {
	skip_on_cran()
	sim = sim_srm(n_actors = 8, seed = 6886)
	pt = permute_srm(sim$Y, n_perms = 30, seed = 6886)
	out = capture.output(print(pt))
	expect_true(any(grepl("Permutation Test", out)))
})

test_that("permute_srm plot returns ggplot", {
	skip_on_cran()
	sim = sim_srm(n_actors = 8, seed = 6886)
	pt = permute_srm(sim$Y, n_perms = 30, seed = 6886)
	p = plot(pt)
	expect_true(inherits(p, "ggplot"))
})

test_that("permute_srm with strong signal has low p-values", {
	skip_on_cran()
	sim = sim_srm(n_actors = 30, actor_var = 8, partner_var = 8,
	              unique_var = 0.3, seed = 6886)
	pt = permute_srm(sim$Y, n_perms = 500, seed = 6886)

	# with very strong actor/partner variance relative to unique,
	# those p-values should be well below 0.05
	expect_true(pt$p_values["actor_var"] < 0.05)
	expect_true(pt$p_values["partner_var"] < 0.05)
})

test_that("permute_srm handles small (n=3) matrices gracefully", {
	skip_on_cran()
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	pt = permute_srm(mat, n_perms = 30, seed = 6886)
	expect_s3_class(pt, "srm_permtest")
	expect_true(all(is.na(pt$p_values)))
})

test_that("permute_srm rejects bipartite srm objects", {
	sim = sim_srm(n_actors = c(6, 8), bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)

	expect_error(permute_srm(fit, n_perms = 5), "unipartite")
})

test_that("permute_srm rejects rectangular matrices", {
	set.seed(6886)
	mat = matrix(rnorm(24), 4, 6)
	rownames(mat) = paste0("S", 1:4)
	colnames(mat) = paste0("R", 1:6)

	expect_error(permute_srm(mat, n_perms = 5), "unipartite")
})

test_that("permute_srm is reproducible with same seed", {
	skip_on_cran()
	sim = sim_srm(n_actors = 8, seed = 6886)
	pt1 = permute_srm(sim$Y, n_perms = 50, seed = 6886)
	pt2 = permute_srm(sim$Y, n_perms = 50, seed = 6886)

	expect_equal(pt1$perm_dist, pt2$perm_dist)
	expect_equal(pt1$p_values, pt2$p_values)
})

test_that("permute_srm works with longitudinal input", {
	skip_on_cran()
	sim = sim_srm(n_actors = 8, n_time = 3, seed = 6886)
	pt = permute_srm(sim$Y, n_perms = 50, seed = 6886)

	expect_s3_class(pt, "srm_permtest")
	expect_equal(pt$n_perms, 50)
	expect_equal(length(pt$observed), 5)
	expect_true(all(pt$p_values >= 0 & pt$p_values <= 1))
})

test_that("permute_srm time filter works on srm objects", {
	skip_on_cran()
	sim = sim_srm(n_actors = 8, n_time = 3, seed = 6886)
	fit = srm(sim$Y)

	# filter via raw matrices
	pt_raw = permute_srm(sim$Y, n_perms = 50, seed = 6886, time = "t2")

	# filter via srm object should give same result
	pt_fit = permute_srm(fit, n_perms = 50, seed = 6886, time = "t2")

	expect_equal(pt_raw$observed, pt_fit$observed, tolerance = 1e-10)
	expect_equal(pt_raw$p_values, pt_fit$p_values)
})

test_that("permute_srm rejects invalid time on srm objects", {
	sim = sim_srm(n_actors = 8, n_time = 2, seed = 6886)
	fit = srm(sim$Y)

	expect_error(permute_srm(fit, n_perms = 5, time = "t99"), "not found")
})

test_that("plot.srm_permtest works with n=3 all-NA results", {
	skip_on_cran()
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")
	pt = permute_srm(mat, n_perms = 20, seed = 6886)

	expect_true(all(is.na(pt$p_values)))
	p = plot(pt)
	expect_true(inherits(p, "ggplot"))
})

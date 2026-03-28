test_that("srm() handles 3x3 minimum size", {
	mat = matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	fit = srm(mat)
	expect_s3_class(fit, "srm")
})

test_that("srm() returns NA variance components for n=3", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	fit = srm(mat)
	expect_true(all(is.na(fit$stats$variance)))
})

test_that("srm() returns valid variance components for n=4", {
	mat = matrix(c(0, 3, 5, 2,
	               7, 0, 4, 6,
	               1, 8, 0, 3,
	               4, 2, 7, 0), 4, 4, byrow = TRUE)
	rownames(mat) = colnames(mat) = c("A", "B", "C", "D")

	fit = srm(mat)
	expect_true(all(is.finite(fit$stats$variance)))
})

test_that("srm() handles zero matrix", {
	mat = matrix(0, 5, 5)
	rownames(mat) = colnames(mat) = letters[1:5]

	fit = srm(mat)
	expect_equal(unname(fit$grand_mean[1]), 0)
	expect_true(all(fit$actor_effects[[1]] == 0))
	expect_true(all(fit$partner_effects[[1]] == 0))
})

test_that("srm() handles large values", {
	mat = matrix(c(0, 1000, 2000, 500, 0, 1500, 3000, 100, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	fit = srm(mat)
	expect_s3_class(fit, "srm")
	expect_true(is.finite(fit$grand_mean[1]))
})

test_that("srm() handles negative values", {
	mat = matrix(c(0, -1, 2, -2, 0, 1, 1, -2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	fit = srm(mat)
	expect_s3_class(fit, "srm")
})

test_that("srm() sets diagonal to zero", {
	mat = matrix(c(5, 1, 2, 2, 5, 1, 1, 2, 5), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	fit = srm(mat)
	expect_equal(unname(diag(fit$unique_effects[[1]])), c(0, 0, 0))

	# grand mean should ignore diagonal values
	mat0 = mat
	diag(mat0) = 0
	fit0 = srm(mat0)
	expect_equal(fit$grand_mean, fit0$grand_mean)
})

test_that("time filtering rejects invalid times", {
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(m1) = colnames(m1) = c("A", "B", "C")

	expect_error(
		srm(list("2020" = m1), time = "1999"),
		"not found"
	)
})

test_that("time filtering on unnamed list fails gracefully", {
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(m1) = colnames(m1) = c("A", "B", "C")

	expect_error(
		srm(list(m1), time = "2020"),
		"no names"
	)
})

test_that("srm() handles matrix without row/colnames", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)

	fit = srm(mat)
	expect_s3_class(fit, "srm")
})

test_that("srm() handles 2x2 matrix", {
	mat = matrix(c(0, 3, 5, 0), 2, 2)
	rownames(mat) = colnames(mat) = c("A", "B")

	fit = srm(mat)
	expect_s3_class(fit, "srm")
	# effects should be NA for n < 3
	expect_true(all(is.na(fit$actor_effects[[1]])))
	expect_true(all(is.na(fit$partner_effects[[1]])))
})

test_that("srm() handles NA values in matrix", {
	mat_na = matrix(c(0, 1, NA, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat_na) = colnames(mat_na) = c("A", "B", "C")

	# NAs are replaced with zero internally
	fit_na = srm(mat_na)
	expect_s3_class(fit_na, "srm")
	expect_true(all(is.finite(fit_na$grand_mean)))

	# result should match a matrix with 0 in place of NA
	mat_zero = mat_na
	mat_zero[is.na(mat_zero)] = 0
	fit_zero = srm(mat_zero)
	expect_equal(fit_na$grand_mean, fit_zero$grand_mean)
	expect_equal(fit_na$actor_effects[[1]], fit_zero$actor_effects[[1]])
})

test_that("srm_effects() assigns default names to unnamed matrices", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	ae = srm_effects(mat, type = "actor")
	expect_equal(rownames(ae), c("n01", "n02", "n03"))

	ue = srm_effects(mat, type = "unique")
	expect_equal(rownames(ue), c("n01", "n02", "n03"))
	expect_equal(colnames(ue), c("n01", "n02", "n03"))
})

test_that("srm_stats() assigns default names to unnamed matrices", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rm = srm_stats(mat, type = "rowmeans")
	expect_equal(names(rm), c("n01", "n02", "n03"))
})

test_that("srm_plot() works on effects from unnamed matrices", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	ae = srm_effects(mat, type = "actor")
	p = srm_plot(ae, type = "actor")
	expect_true(inherits(p, "ggplot"))
	expect_equal(nrow(p$data), 3)
})

test_that("srm_stats() handles 1x1 matrix without NaN", {
	mat = matrix(0, 1, 1)
	rownames(mat) = colnames(mat) = "A"

	rm = srm_stats(mat, type = "rowmeans")
	expect_true(is.na(rm))

	tm = srm_stats(mat, type = "totalmeans")
	expect_true(is.na(tm))

	av = srm_stats(mat, type = "actor_var")
	expect_true(is.na(av))
})

test_that("srm_effects() gives informative error for rectangular matrix", {
	mat = matrix(1:6, 2, 3)
	expect_error(srm_effects(mat, type = "actor"), "bipartite")
})

test_that("srm_stats() gives informative error for rectangular matrix", {
	mat = matrix(1:6, 2, 3)
	expect_error(srm_stats(mat, type = "rowmeans"), "bipartite")
})

test_that("sim_srm() handles zero variances", {
	sim = sim_srm(n_actors = 5, actor_var = 0, partner_var = 0,
	              unique_var = 0, seed = 6886)
	expect_equal(unname(diag(sim$Y)), rep(0, 5))
	off_diag = sim$Y[row(sim$Y) != col(sim$Y)]
	expect_true(all(off_diag == 0))

	# with nonzero grand mean
	sim2 = sim_srm(n_actors = 5, actor_var = 0, partner_var = 0,
	               unique_var = 0, grand_mean = 3.0, seed = 6886)
	off_diag2 = sim2$Y[row(sim2$Y) != col(sim2$Y)]
	expect_true(all(off_diag2 == 3.0))
})

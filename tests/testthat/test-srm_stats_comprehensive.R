test_that("srm_stats rowmeans is correct", {
	# matrix(c(0,2,4, 6,0,8, 10,12,0), 3, 3) fills by column:
	# col1: 0,2,4 | col2: 6,0,8 | col3: 10,12,0
	# row a: 0, 6, 10 -> mean = (0+6+10)/2 = 8
	mat = matrix(c(0, 2, 4, 6, 0, 8, 10, 12, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	rm = srm_stats(mat, type = "rowmeans")
	expect_equal(length(rm), 3)
	# after diag set to 0: row a = (0, 6, 10), mean = 16/2 = 8
	expect_equal(rm["A"], c(A = (6 + 10) / 2))
})

test_that("srm_stats colmeans is correct", {
	# matrix fills by column: col1 = (0,2,4), col2 = (6,0,8), col3 = (10,12,0)
	# after diag set to 0: col A = (0,2,4), mean = (2+4)/2 = 3
	# col B = (6,0,8), mean = (6+8)/2 = 7
	# col C = (10,12,0), mean = (10+12)/2 = 11
	mat = matrix(c(0, 2, 4, 6, 0, 8, 10, 12, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	cm = srm_stats(mat, type = "colmeans")
	expect_equal(length(cm), 3)
	expect_equal(cm["A"], c(A = (2 + 4) / 2))
	expect_equal(cm["B"], c(B = (6 + 8) / 2))
	expect_equal(cm["C"], c(C = (10 + 12) / 2))
})

test_that("srm_stats totalmeans is correct", {
	mat = matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	tm = srm_stats(mat, type = "totalmeans")
	expect_equal(tm, sum(mat) / (3 * 2))
})

test_that("srm_stats all variance types run without error", {
	mat = matrix(c(0, 3, 5, 7, 0, 2, 4, 6, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	types = c("actor_var", "partner_var", "unique_var",
	          "relationship_cov", "actor_partner_cov")

	for (tp in types) {
		result = srm_stats(mat, type = tp)
		expect_true(is.numeric(result), info = paste("Failed for type:", tp))
		expect_equal(length(result), 1, info = paste("Failed for type:", tp))
	}
})

test_that("srm_stats longitudinal returns list", {
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = m1 * 2
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	rm = srm_stats(list("t1" = m1, "t2" = m2), type = "rowmeans")
	expect_true(is.list(rm))
	expect_equal(length(rm), 2)
})

test_that("srm_stats rejects non-square matrix", {
	mat = matrix(1:6, 2, 3)
	expect_error(srm_stats(mat, type = "rowmeans"), "square")
})

test_that("srm_stats variance components agree with srm()", {
	sim = sim_srm(n_actors = 8, seed = 6886)
	mat = sim$Y

	fit = srm(mat)
	av_srm = fit$stats$variance[fit$stats$component == "actor_var"]
	av_stats = srm_stats(mat, type = "actor_var")

	expect_equal(av_srm, av_stats, tolerance = 1e-10)
})

test_that("srm_stats returns NA for variance types with n=3", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	types = c("actor_var", "partner_var", "unique_var",
	          "relationship_cov", "actor_partner_cov")

	for (tp in types) {
		result = srm_stats(mat, type = tp)
		expect_true(is.na(result), info = paste("Expected NA for type:", tp, "with n=3"))
	}
})

test_that("srm_stats all five components agree with srm() for all types", {
	sim = sim_srm(n_actors = 10, seed = 6886)
	mat = sim$Y
	fit = srm(mat)

	types = c("actor_var", "partner_var", "unique_var",
	          "relationship_cov", "actor_partner_cov")

	for (tp in types) {
		from_srm = fit$stats$variance[fit$stats$component == tp]
		from_stats = srm_stats(mat, type = tp)
		expect_equal(from_srm, from_stats, tolerance = 1e-10,
		             info = paste("Mismatch for:", tp))
	}
})

test_that("srm_stats rejects srm objects with helpful message", {
	sim = sim_srm(n_actors = 6, seed = 6886)
	fit = srm(sim$Y)

	expect_error(srm_stats(fit, type = "rowmeans"), "srm")
})

test_that("srm_stats handles netify objects", {
	skip_if_not_installed("netify")
	data(classroom)
	# compute from raw matrix
	rm_raw = srm_stats(classroom, type = "rowmeans")

	# wrap as netify-like matrix and compute
	mat_net = classroom
	class(mat_net) = c("netify", "matrix")
	rm_net = suppressMessages(srm_stats(mat_net, type = "rowmeans"))

	expect_equal(rm_raw, rm_net)
})

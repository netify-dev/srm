test_that("srm_trends returns tidy data frame with correct values", {
	sim = sim_srm(n_actors = 5, n_time = 3, seed = 6886)
	fit = srm(sim$Y)

	trends = srm_trends(fit, type = "actor")

	expect_true(is.data.frame(trends))
	expect_true(all(c("actor", "time", "effect") %in% names(trends)))
	expect_equal(nrow(trends), 5 * 3)

	# values should match actor effects stored in the srm object
	for (tp in names(fit$actor_effects)) {
		eff = fit$actor_effects[[tp]]
		sub = trends[trends$time == tp, ]
		for (a in rownames(eff)) {
			expect_equal(
				sub$effect[sub$actor == a],
				eff[a, 1],
				tolerance = 1e-10,
				info = paste("mismatch for", a, "at", tp)
			)
		}
	}
})

test_that("srm_trends filters actors", {
	sim = sim_srm(n_actors = 6, n_time = 3, seed = 6886)
	fit = srm(sim$Y)

	trends = srm_trends(fit, type = "actor", actors = c("n01", "n02"))
	expect_true(all(trends$actor %in% c("n01", "n02")))
})

test_that("srm_trends works for partner type", {
	sim = sim_srm(n_actors = 5, n_time = 2, seed = 6886)
	fit = srm(sim$Y)

	trends = srm_trends(fit, type = "partner")
	expect_equal(nrow(trends), 5 * 2)
})

test_that("srm_trend_plot returns ggplot", {
	skip_on_cran()
	sim = sim_srm(n_actors = 6, n_time = 3, seed = 6886)
	fit = srm(sim$Y)

	p = srm_trend_plot(fit, type = "actor", n = 3)
	expect_true(inherits(p, "ggplot"))
})

test_that("srm_stability returns correct correlations", {
	skip_on_cran()
	sim = sim_srm(n_actors = 10, n_time = 4, seed = 6886)
	fit = srm(sim$Y)

	stab = srm_stability(fit, type = "actor")

	expect_true(is.data.frame(stab))
	expect_equal(nrow(stab), 3)
	expect_true(all(c("time1", "time2", "correlation", "n") %in% names(stab)))
	expect_true(all(stab$correlation >= -1 & stab$correlation <= 1))

	# verify first correlation by hand
	e1 = fit$actor_effects[["t1"]][, 1]
	e2 = fit$actor_effects[["t2"]][, 1]
	shared = intersect(names(e1), names(e2))
	expected_r = cor(e1[shared], e2[shared])
	expect_equal(stab$correlation[1], expected_r, tolerance = 1e-10)
	expect_equal(stab$n[1], length(shared))
})

test_that("srm_stability works for partner type", {
	skip_on_cran()
	sim = sim_srm(n_actors = 10, n_time = 4, seed = 6886)
	fit = srm(sim$Y)

	stab = srm_stability(fit, type = "partner")

	expect_true(is.data.frame(stab))
	expect_equal(nrow(stab), 3)
	expect_true(all(stab$correlation >= -1 & stab$correlation <= 1))
})

test_that("srm_stability errors on single time point", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")
	fit = srm(mat)

	expect_error(srm_stability(fit), "at least 2")
})

test_that("srm_trends errors on non-srm input", {
	expect_error(srm_trends("not srm"), class = "rlang_error")
})

test_that("srm_trend_plot works with specific actors", {
	skip_on_cran()
	sim = sim_srm(n_actors = 6, n_time = 3, seed = 6886)
	fit = srm(sim$Y)

	p = srm_trend_plot(fit, type = "actor", actors = c("n01", "n02"))
	expect_true(inherits(p, "ggplot"))
})

test_that("plot.srm time filtering works for longitudinal", {
	skip_on_cran()
	sim = sim_srm(n_actors = 6, n_time = 3, seed = 6886)
	fit = srm(sim$Y)

	p = plot(fit, type = "actor", time = "t2")
	expect_true(inherits(p, "ggplot"))

	p_var = plot(fit, type = "variance", time = c("t1", "t3"))
	expect_true(inherits(p_var, "ggplot"))
})

test_that("srm_trends works for bipartite longitudinal", {
	skip_on_cran()
	sim = sim_srm(n_actors = c(5, 7), n_time = 3, bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)

	trends_a = srm_trends(fit, type = "actor")
	expect_equal(nrow(trends_a), 5 * 3)

	trends_p = srm_trends(fit, type = "partner")
	expect_equal(nrow(trends_p), 7 * 3)
})

test_that("srm_stability works for bipartite longitudinal", {
	skip_on_cran()
	sim = sim_srm(n_actors = c(5, 7), n_time = 3, bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)

	stab_a = srm_stability(fit, type = "actor")
	expect_equal(nrow(stab_a), 2)
	expect_true(all(stab_a$n == 5))

	stab_p = srm_stability(fit, type = "partner")
	expect_equal(nrow(stab_p), 2)
	expect_true(all(stab_p$n == 7))
})

test_that("srm_trend_plot works for bipartite longitudinal", {
	skip_on_cran()
	sim = sim_srm(n_actors = c(5, 7), n_time = 3, bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)

	p = srm_trend_plot(fit, type = "actor", n = 3)
	expect_true(inherits(p, "ggplot"))

	p2 = srm_trend_plot(fit, type = "partner", n = 3)
	expect_true(inherits(p2, "ggplot"))
})

test_that("plot.srm time filtering works for partner and dyadic types", {
	skip_on_cran()
	sim = sim_srm(n_actors = 6, n_time = 3, seed = 6886)
	fit = srm(sim$Y)

	p_partner = plot(fit, type = "partner", time = "t2")
	expect_true(inherits(p_partner, "ggplot"))

	p_dyadic = plot(fit, type = "dyadic", time = "t1")
	expect_true(inherits(p_dyadic, "ggplot"))
})

test_that("srm_stability returns NA when fewer than 3 shared actors", {
	# two periods with only 2 overlapping actors
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = matrix(c(0, 3, 1, 1, 0, 2, 2, 1, 0), 3, 3)
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "D")

	fit = srm(list("t1" = m1, "t2" = m2))
	stab = srm_stability(fit, type = "actor")

	expect_equal(nrow(stab), 1)
	expect_true(is.na(stab$correlation[1]))
	expect_equal(stab$n[1], 2)
})

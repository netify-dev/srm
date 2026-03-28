test_that("print.srm works", {
	sim = sim_srm(n_actors = 8, seed = 6886)
	fit = srm(sim$Y)

	out = capture.output(print(fit))
	expect_true(any(grepl("Social Relations Model", out)))
	expect_true(any(grepl("Actor", out)))
	expect_true(any(grepl("Partner", out)))
	expect_true(any(grepl("Unique", out)))
})

test_that("print.srm works on small (n=3) matrix", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")
	fit = srm(mat)

	out = capture.output(print(fit))
	expect_true(any(grepl("Social Relations Model", out)))
	expect_true(any(grepl("too small", out)))
})

test_that("summary.srm returns correct class", {
	sim = sim_srm(n_actors = 8, seed = 6886)
	fit = srm(sim$Y)
	s = summary(fit)

	expect_s3_class(s, "summary.srm")
	expect_true("variance_table" %in% names(s))
	expect_equal(nrow(s$variance_table), 5) # 5 components
})

test_that("print.summary.srm works", {
	sim = sim_srm(n_actors = 8, seed = 6886)
	fit = srm(sim$Y)
	s = summary(fit)

	out = capture.output(print(s))
	expect_true(any(grepl("Variance Decomposition", out)))
	expect_true(any(grepl("Actor", out)))
})

test_that("summary variance percentages sum to 100", {
	sim = sim_srm(n_actors = 12, seed = 6886)
	fit = srm(sim$Y)
	s = summary(fit)

	var_pct = s$variance_table$pct[!is.na(s$variance_table$pct)]
	expect_equal(length(var_pct), 3)
	expect_equal(sum(var_pct), 100, tolerance = 0.01)
})

test_that("bipartite summary percentages sum to 100", {
	sim = sim_srm(n_actors = c(8, 10), bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)
	s = summary(fit)

	var_pct = s$variance_table$pct[!is.na(s$variance_table$pct)]
	expect_equal(length(var_pct), 3)
	expect_equal(sum(var_pct), 100, tolerance = 0.01)
})

test_that("plot.srm returns ggplot objects", {
	sim = sim_srm(n_actors = 8, seed = 6886)
	fit = srm(sim$Y)

	p_actor = plot(fit, type = "actor")
	expect_true(inherits(p_actor, "ggplot"))

	p_partner = plot(fit, type = "partner")
	expect_true(inherits(p_partner, "ggplot"))

	p_dyadic = plot(fit, type = "dyadic")
	expect_true(inherits(p_dyadic, "ggplot"))

	p_var = plot(fit, type = "variance")
	expect_true(inherits(p_var, "ggplot"))
})

test_that("print.srm shows longitudinal info", {
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = m1 * 2
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	fit = srm(list("2020" = m1, "2021" = m2))
	out = capture.output(print(fit))
	expect_true(any(grepl("Time points", out)))
})

test_that("negative variance shows -- instead of 0.0% in print", {
	# matrix that produces negative actor/partner variance
	mat = matrix(c(0, 10, 1, 2,
	               1, 0, 10, 2,
	               2, 1, 0, 10,
	               10, 2, 1, 0), 4, 4, byrow = TRUE)
	rownames(mat) = colnames(mat) = LETTERS[1:4]
	fit = srm(mat)

	# confirm negative variance exists
	expect_true(any(fit$stats$variance < 0))

	# print should show -- for negative variance percentages
	out = capture.output(print(fit))
	actor_line = out[grepl("Actor", out) & grepl("--", out)]
	expect_true(length(actor_line) > 0)

	# summary should also show -- for negative variance
	out_s = capture.output(print(summary(fit)))
	actor_line_s = out_s[grepl("Actor", out_s) & !grepl("Partner", out_s)]
	expect_true(any(grepl("--", actor_line_s)))
})

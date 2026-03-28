test_that("srm_plot actor works with small matrix", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")
	ae = srm_effects(mat, type = "actor")

	p = srm_plot(ae, type = "actor", n = 3)
	expect_true(inherits(p, "ggplot"))
})

test_that("srm_plot partner works", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")
	pe = srm_effects(mat, type = "partner")

	p = srm_plot(pe, type = "partner", n = 3)
	expect_true(inherits(p, "ggplot"))
})

test_that("srm_plot dyadic works", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")
	ue = srm_effects(mat, type = "unique")

	p = srm_plot(ue, type = "dyadic", n = 3)
	expect_true(inherits(p, "ggplot"))
})

test_that("srm_plot facet works for longitudinal actor", {
	skip_on_cran()
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = m1 * 2
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	ae = srm_effects(list("t1" = m1, "t2" = m2), type = "actor")
	p = srm_plot(ae, type = "actor", facet = TRUE, n = 3)
	expect_true(inherits(p, "ggplot"))
})

test_that("srm_plot dyadic facet gives warning", {
	skip_on_cran()
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = m1 * 2
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	ue = srm_effects(list("t1" = m1, "t2" = m2), type = "unique")
	expect_warning(srm_plot(ue, type = "dyadic", facet = TRUE, n = 3))
})

test_that("srm_plot n parameter limits output", {
	skip_on_cran()
	set.seed(6886)
	n = 10
	mat = matrix(rpois(n * n, 3), n, n)
	diag(mat) = 0
	rownames(mat) = colnames(mat) = paste0("n", sprintf("%02d", 1:n))

	ae = srm_effects(mat, type = "actor")
	p = srm_plot(ae, type = "actor", n = 5)
	expect_true(inherits(p, "ggplot"))

	# the plot data should have at most 5 rows
	expect_true(nrow(p$data) <= 5)
})

test_that("srm_plot dyadic works for bipartite", {
	skip_on_cran()
	sim = sim_srm(n_actors = c(6, 8), bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)

	p = plot(fit, type = "dyadic", n = 5)
	expect_true(inherits(p, "ggplot"))
})

test_that("srm_plot actor works for bipartite", {
	skip_on_cran()
	sim = sim_srm(n_actors = c(6, 8), bipartite = TRUE, seed = 6886)
	fit = srm(sim$Y)

	p = plot(fit, type = "actor")
	expect_true(inherits(p, "ggplot"))
})

test_that("srm_plot rejects bad input", {
	expect_error(srm_plot("not a matrix", type = "actor"))
})

test_that("srm_plot time filter works", {
	skip_on_cran()
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = m1 * 2
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	ae = srm_effects(list("2020" = m1, "2021" = m2), type = "actor")
	p = srm_plot(ae, type = "actor", time = "2020", n = 3)
	expect_true(inherits(p, "ggplot"))
})

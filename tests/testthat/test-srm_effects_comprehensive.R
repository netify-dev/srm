test_that("srm_effects actor returns correct structure", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	ae = srm_effects(mat, type = "actor")
	expect_equal(nrow(ae), 3)
	expect_equal(ncol(ae), 1)
	expect_equal(rownames(ae), c("A", "B", "C"))
	expect_true(inherits(ae, "actor_effect"))
})

test_that("srm_effects partner returns correct structure", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	pe = srm_effects(mat, type = "partner")
	expect_equal(nrow(pe), 3)
	expect_equal(ncol(pe), 1)
	expect_true(inherits(pe, "partner_effect"))
})

test_that("srm_effects unique returns square matrix", {
	mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	rownames(mat) = colnames(mat) = c("A", "B", "C")

	ue = srm_effects(mat, type = "unique")
	expect_equal(nrow(ue), 3)
	expect_equal(ncol(ue), 3)
	expect_true(inherits(ue, "unique_effect"))
	expect_equal(unname(diag(ue)), c(0, 0, 0))
})

test_that("srm_effects works on list of matrices", {
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = m1 * 2
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	ae = srm_effects(list("t1" = m1, "t2" = m2), type = "actor")
	expect_true(is.list(ae))
	expect_equal(length(ae), 2)
})

test_that("srm_effects rejects non-square matrix", {
	mat = matrix(1:6, 2, 3)
	expect_error(srm_effects(mat, type = "actor"), "square")
})

test_that("srm_effects actor sum is approximately zero", {
	set.seed(6886)
	n = 10
	mat = matrix(rpois(n * n, 3), n, n)
	diag(mat) = 0
	rownames(mat) = colnames(mat) = paste0("n", 1:n)

	ae = srm_effects(mat, type = "actor")
	# actor effects should approximately sum to zero
	expect_true(abs(sum(ae[, 1])) < 1e-10)
})

test_that("srm_effects time filtering works", {
	m1 = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
	m2 = m1 * 3
	rownames(m1) = colnames(m1) = c("A", "B", "C")
	rownames(m2) = colnames(m2) = c("A", "B", "C")

	ae = srm_effects(list("2020" = m1, "2021" = m2), type = "actor", time = "2021")
	expect_true(is.matrix(ae))  # single matrix returned (not list)
})

test_that("srm_effects agrees with srm() decomposition", {
	sim = sim_srm(n_actors = 8, seed = 6886)
	mat = sim$Y
	fit = srm(mat)

	ae = srm_effects(mat, type = "actor")
	pe = srm_effects(mat, type = "partner")
	ue = srm_effects(mat, type = "unique")

	expect_equal(ae[, 1], fit$actor_effects[[1]][, 1], tolerance = 1e-10)
	expect_equal(pe[, 1], fit$partner_effects[[1]][, 1], tolerance = 1e-10)
	expect_equal(unname(ue), unname(fit$unique_effects[[1]]), tolerance = 1e-10)
})

test_that("srm_effects decomposition reconstructs original off-diagonals", {
	sim = sim_srm(n_actors = 8, seed = 6886)
	mat = sim$Y
	diag(mat) = 0
	n = nrow(mat)

	ae = srm_effects(mat, type = "actor")
	pe = srm_effects(mat, type = "partner")
	ue = srm_effects(mat, type = "unique")
	mu = sum(mat) / (n * (n - 1))

	reconstructed = mu +
		matrix(ae[, 1], n, n) +
		matrix(pe[, 1], n, n, byrow = TRUE) +
		ue

	for (i in 1:n) {
		for (j in 1:n) {
			if (i != j) {
				expect_equal(reconstructed[i, j], mat[i, j], tolerance = 1e-10)
			}
		}
	}
})

test_that("srm_effects column names match srm() output", {
	sim = sim_srm(n_actors = 6, seed = 6886)
	mat = sim$Y
	fit = srm(mat)

	ae = srm_effects(mat, type = "actor")
	pe = srm_effects(mat, type = "partner")

	expect_equal(colnames(ae), colnames(fit$actor_effects[[1]]))
	expect_equal(colnames(pe), colnames(fit$partner_effects[[1]]))
})

test_that("srm_effects rejects srm objects with helpful message", {
	sim = sim_srm(n_actors = 6, seed = 6886)
	fit = srm(sim$Y)

	expect_error(srm_effects(fit, type = "actor"), "srm")
})

test_that("srm_effects returns NA for matrices with fewer than 3 actors", {
	mat = matrix(c(0, 1, 2, 0), 2, 2)
	rownames(mat) = colnames(mat) = c("A", "B")

	ae = srm_effects(mat, type = "actor")
	expect_true(all(is.na(ae[, 1])))
	expect_true(inherits(ae, "actor_effect"))
	expect_equal(nrow(ae), 2)

	pe = srm_effects(mat, type = "partner")
	expect_true(all(is.na(pe[, 1])))
	expect_true(inherits(pe, "partner_effect"))

	ue = srm_effects(mat, type = "unique")
	expect_true(inherits(ue, "unique_effect"))
	expect_equal(diag(ue), c(A = 0, B = 0))
	expect_true(all(is.na(ue[row(ue) != col(ue)])))
})

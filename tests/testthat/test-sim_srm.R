test_that("sim_srm generates correct dimensions", {
	sim = sim_srm(n_actors = 8, seed = 6886)

	expect_true(is.matrix(sim$Y))
	expect_equal(nrow(sim$Y), 8)
	expect_equal(ncol(sim$Y), 8)
	expect_equal(length(sim$truth$actor_effects), 8)
	expect_equal(length(sim$truth$partner_effects), 8)

	# diagonal should be zero
	expect_equal(unname(diag(sim$Y)), rep(0, 8))
})

test_that("sim_srm longitudinal generates list", {
	sim = sim_srm(n_actors = 6, n_time = 3, seed = 6886)

	expect_true(is.list(sim$Y))
	expect_equal(length(sim$Y), 3)
	expect_equal(names(sim$Y), c("t1", "t2", "t3"))
	expect_equal(nrow(sim$Y$t1), 6)
})

test_that("sim_srm bipartite works", {
	sim = sim_srm(n_actors = c(5, 8), bipartite = TRUE, seed = 6886)

	expect_equal(nrow(sim$Y), 5)
	expect_equal(ncol(sim$Y), 8)
	expect_true(sim$params$bipartite)
})

test_that("sim_srm bipartite longitudinal works", {
	sim = sim_srm(n_actors = c(4, 6), n_time = 3, bipartite = TRUE, seed = 6886)

	expect_true(is.list(sim$Y))
	expect_equal(length(sim$Y), 3)
	expect_equal(nrow(sim$Y$t1), 4)
	expect_equal(ncol(sim$Y$t1), 6)
})

test_that("sim_srm seed is reproducible", {
	sim1 = sim_srm(n_actors = 5, seed = 6886)
	sim2 = sim_srm(n_actors = 5, seed = 6886)
	expect_equal(sim1$Y, sim2$Y)
})

test_that("sim_srm with srm() recovers approximate variances", {
	skip_on_cran()
	sim = sim_srm(n_actors = 30, actor_var = 2, partner_var = 1,
	              unique_var = 3, seed = 6886)
	fit = srm(sim$Y)
	s = summary(fit)

	# large n should give reasonable estimates (within factor of 2)
	vt = s$variance_table
	est_actor = vt$variance[vt$component == "actor_var"]
	est_partner = vt$variance[vt$component == "partner_var"]
	est_unique = vt$variance[vt$component == "unique_var"]

	expect_true(est_actor > 0.5 && est_actor < 6)
	expect_true(est_partner > 0.2 && est_partner < 4)
	expect_true(est_unique > 1 && est_unique < 8)
})

test_that("sim_srm bipartite recovers approximate variances", {
	skip_on_cran()
	sim = sim_srm(n_actors = c(20, 25), bipartite = TRUE,
	              actor_var = 2, partner_var = 1, unique_var = 1.5,
	              seed = 6886)
	fit = srm(sim$Y)
	s = summary(fit)

	vt = s$variance_table
	est_actor = vt$variance[vt$component == "actor_var"]
	est_partner = vt$variance[vt$component == "partner_var"]
	est_unique = vt$variance[vt$component == "unique_var"]

	# reasonable recovery for bipartite with moderate n
	expect_true(est_actor > 0.5 && est_actor < 6)
	expect_true(est_partner > 0.2 && est_partner < 4)
	expect_true(est_unique > 0.4 && est_unique < 5)
})

test_that("sim_srm rejects invalid n_actors", {
	expect_error(sim_srm(n_actors = 2), "n_actors")
})

test_that("sim_srm rejects invalid actor-partner covariance", {
	expect_error(sim_srm(n_actors = 5, actor_var = 1, partner_var = 1,
	                     actor_partner_cov = 5),
	             "positive semi-definite")
})

test_that("sim_srm rejects invalid relationship covariance", {
	# relationship_cov > unique_var is not PSD
	expect_error(sim_srm(n_actors = 5, unique_var = 1, relationship_cov = 2),
	             "positive semi-definite")
	# negative overflow
	expect_error(sim_srm(n_actors = 5, unique_var = 1, relationship_cov = -2),
	             "positive semi-definite")
})

test_that("sim_srm handles boundary relationship covariance", {
	# relationship_cov = unique_var is PSD (perfect reciprocity)
	sim = sim_srm(n_actors = 10, unique_var = 1, relationship_cov = 1, seed = 6886)
	expect_true(is.matrix(sim$Y))
	expect_equal(nrow(sim$Y), 10)
})

test_that("sim_srm params are stored", {
	sim = sim_srm(n_actors = 5, actor_var = 2.5, seed = 6886)
	expect_equal(sim$params$actor_var, 2.5)
	expect_equal(sim$params$n_actors, 5)
	expect_false(sim$params$bipartite)
})

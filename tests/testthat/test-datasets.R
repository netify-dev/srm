test_that("atop_EA dataset loads correctly", {
	data(atop_EA, package = "srm")
	expect_true(is.data.frame(atop_EA))
	expect_true(all(c("year", "country1", "country2") %in% names(atop_EA)))
	expect_true(nrow(atop_EA) > 0)
})

test_that("classroom dataset loads and fits", {
	data(classroom, package = "srm")

	expect_true(is.matrix(classroom))
	expect_equal(nrow(classroom), 12)
	expect_equal(ncol(classroom), 12)
	expect_equal(unname(diag(classroom)), rep(0, 12))

	# mean girls character names
	expect_true("Regina" %in% rownames(classroom))
	expect_true("Cady" %in% rownames(classroom))

	fit = srm(classroom)
	expect_s3_class(fit, "srm")
})

test_that("classroom decomposition matches vignette values", {
	# regression test: vignette narratives reference these exact values
	data(classroom, package = "srm")
	fit = srm(classroom)

	# actor effects used in pipeline vignette
	a = fit$actor_effects[[1]][, 1]
	expect_equal(unname(round(a["Damian"], 2)), 2.36)
	expect_equal(unname(round(a["Regina"], 2)), -2.64)

	# partner effects used in pipeline vignette
	b = fit$partner_effects[[1]][, 1]
	expect_equal(unname(round(b["Regina"], 2)), 1.67)
	expect_equal(unname(round(b["Gretchen"], 2)), -1.71)

	# variance decomposition percentages
	s = summary(fit)
	vt = s$variance_table
	var_pct = vt$pct[!is.na(vt$pct)]
	expect_equal(round(sum(var_pct), 1), 100)
	expect_equal(round(var_pct[1], 1), 42.3)  # actor
	expect_equal(round(var_pct[2], 1), 25.8)  # partner
	expect_equal(round(var_pct[3], 1), 31.9)  # unique
})

test_that("trade_net dataset loads and fits", {
	data(trade_net, package = "srm")

	expect_true(is.list(trade_net))
	expect_equal(length(trade_net), 3)
	expect_equal(names(trade_net), c("2015", "2017", "2019"))

	fit = srm(trade_net)
	expect_s3_class(fit, "srm")
	expect_equal(fit$n_time, 3)
})

test_that("small_council dataset loads and fits as bipartite", {
	data(small_council, package = "srm")

	expect_true(is.matrix(small_council))
	expect_equal(nrow(small_council), 10)
	expect_equal(ncol(small_council), 7)
	expect_true("Lannister" %in% rownames(small_council))
	expect_true("Hand" %in% colnames(small_council))

	fit = srm(small_council)
	expect_s3_class(fit, "srm")
	expect_true(fit$bipartite)
})

test_that("small_council decomposition matches vignette values", {
	data(small_council, package = "srm")
	fit = srm(small_council)

	# actor effects referenced in bipartite vignette
	a = fit$actor_effects[[1]][, 1]
	expect_equal(unname(round(a["Lannister"], 2)), 3.09)
	expect_equal(unname(round(a["Arryn"], 2)), -1.88)

	# partner effects referenced in bipartite vignette
	b = fit$partner_effects[[1]][, 1]
	expect_equal(unname(round(b["Hand"], 2)), 1.63)
	expect_equal(unname(round(b["Faith"], 2)), -2.20)

	# unique effects referenced in bipartite vignette
	u = fit$unique_effects[[1]]
	expect_equal(unname(round(u["Greyjoy", "Ships"], 2)), 3.08)
	expect_equal(unname(round(u["Stark", "Whispers"], 2)), -2.59)

	# variance partition
	s = summary(fit)
	vt = s$variance_table
	var_pct = vt$pct[!is.na(vt$pct)]
	expect_equal(round(sum(var_pct), 1), 100)
	expect_equal(round(var_pct[1], 1), 37.8)  # actor
	expect_equal(round(var_pct[2], 1), 28.9)  # partner
	expect_equal(round(var_pct[3], 1), 33.2)  # unique
})

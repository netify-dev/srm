test_that("srm_stats works", {
  # Create a sample matrix
  test_matrix <- matrix(c(0, 1, 2, 
                          2, 0, 1, 
                          1, 2, 0), nrow=3, ncol=3)
  rownames(test_matrix) <- c("A", "B", "C")
  colnames(test_matrix) <- c("A", "B", "C")
  
  # Test rowmeans only
  row_means <- srm_stats(test_matrix, type="rowmeans")
  
  # Check output
  expect_equal(length(row_means), 3)
  expect_equal(names(row_means), c("A", "B", "C"))
})
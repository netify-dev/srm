
test_that("srm_effects actor works", {
  # Create a sample matrix
  test_matrix <- matrix(c(0, 1, 2, 
                          2, 0, 1, 
                          1, 2, 0), nrow=3, ncol=3)
  rownames(test_matrix) <- c("A", "B", "C")
  colnames(test_matrix) <- c("A", "B", "C")
  
  # Test actor effects only
  actor_effects <- srm_effects(test_matrix, type="actor")
  
  # Check output dimensions and structure
  expect_equal(nrow(actor_effects), 3)
  expect_equal(ncol(actor_effects), 1)
  expect_equal(rownames(actor_effects), c("A", "B", "C"))
})
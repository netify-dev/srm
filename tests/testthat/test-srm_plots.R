test_that("srm_plot works", {
  # Create a sample matrix
  test_matrix <- matrix(c(0, 1, 2, 
                          2, 0, 1, 
                          1, 2, 0), nrow=3, ncol=3)
  rownames(test_matrix) <- c("A", "B", "C")
  colnames(test_matrix) <- c("A", "B", "C")
  
  # Calculate actor effects
  actor_effects <- srm_effects(test_matrix, type="actor")
  
  # Test actor plot only
  actor_plot <- srm_plot(actor_effects, type="actor", n=3)
  expect_true(inherits(actor_plot, "ggplot"))
})
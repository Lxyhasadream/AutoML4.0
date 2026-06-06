test_that("auto_ml_analysis is exported", {
  expect_true(is.function(auto_ml_analysis))
  expect_true(is.function(aml_methods))
  expect_length(aml_methods(), 25)
})

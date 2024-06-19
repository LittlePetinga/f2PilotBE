library(testthat)
library(f2PilotBE)

# Load data from CSV file
dta <- read.csv(system.file("data", "NormDta.csv", package = "f2PilotBE"))

# Test f2 function
test_that("Test for f2", {

  # Each test validates the computation of f2 factor using different test scenarios (Test1 to Test4).
  # The expected values are based on predefined expectations or empirical data.

  # Test 1: Verify f2 computation for Test1
  Test1 <- f2(dta, Ref = 'Reference', Test = 'Test1')
  expect_equal(round(Test1$f2,0), 52)

  # Test 2: Verify f2 computation for Test2
  Test2 <- f2(dta, Ref = 'Reference', Test = 'Test2')
  expect_equal(round(Test2$f2,0), 72)

  # Test 3: Verify f2 computation for Test3
  Test3 <- f2(dta, Ref = 'Reference', Test = 'Test3')
  expect_equal(round(Test3$f2,0), 59)

  # Test 4: Verify f2 computation for Test4
  Test4 <- f2(dta, Ref = 'Reference', Test = 'Test4')
  expect_equal(round(Test4$f2,0), 36)
})

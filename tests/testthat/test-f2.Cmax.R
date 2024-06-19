library(testthat)
library(f2PilotBE)

# Load data from CSV file
Conc1 <- read.csv(system.file("data", "Conc1.csv", package = "f2PilotBE"))
Conc2 <- read.csv(system.file("data", "Conc2.csv", package = "f2PilotBE"))
Conc3 <- read.csv(system.file("data", "Conc3.csv", package = "f2PilotBE"))
Conc4 <- read.csv(system.file("data", "Conc4.csv", package = "f2PilotBE"))
Conc5 <- read.csv(system.file("data", "Conc5.csv", package = "f2PilotBE"))


# Test f2.Cmax function
test_that("Test for f2.Cmax", {

  # Each test validates the computation of f2 factor using different mean
  # concentration-time profiles of Test and Reference products.
  # The expected values are based on predefined expectations or empirical data.

  # Test 1: Verify Cmax f2 computation for Conc1
  # Using stacked treatment information
  # Where Reference presents an higher bioavailability
  # But within predefined limits
  Test1 <- f2.Cmax(Conc1, Time = 'time', Conc = "conc",
                   Trt = "trt", Ref = 'R', Test = 'T',
                   Trt.cols = FALSE, details = FALSE, plot = FALSE)
  expect_equal(round(Test1$`Cmax f2 Factor`$f2,0), 47)


  # Test 2: Verify Cmax f2 computation for Conc2
  # Using stacked treatment information
  # Where Test and Reference presents a similar bioavailability
  Test2 <- f2.Cmax(Conc2, Time = 'time', Conc = "conc",
                   Trt = "trt", Ref = 'A', Test = 'B',
                   Trt.cols = FALSE, details = FALSE, plot = FALSE)
  expect_equal(round(Test2$`Cmax f2 Factor`$f2,0), 83)


  # Test 3: Verify Cmax f2 computation for Conc3
  # Using stacked treatment information
  # Where Reference presents an higher bioavailability (and not bioequivalent)
  Test3 <- f2.Cmax(Conc3, Time = 'time', Conc = "conc",
                   Trt = "trt", Ref = 'R', Test = 'T',
                   Trt.cols = FALSE, details = FALSE, plot = FALSE)
  expect_equal(round(Test3$`Cmax f2 Factor`$f2,0), 21)


  # Test 4: Verify Cmax f2 computation for Conc4
  # Using stacked treatment information
  # Where Reference presents a lower bioavailability (and not bioequivalent)
  Test4 <- f2.Cmax(Conc4, Time = 'time', Conc = "conc",
                   Trt = "trt", Ref = 'R', Test = 'T',
                   Trt.cols = FALSE, details = FALSE, plot = FALSE)
  expect_equal(round(Test4$`Cmax f2 Factor`$f2,0), 15)


  # Test 5: Verify Cmax f2 computation for Conc5
  # Using stacked treatment information
  # Where Reference presents an higher bioavailability (and not bioequivalent)
  Test5 <- f2.Cmax(Conc5, Time = 'Time',Ref = 'Reference', Test = 'Test',
                   Trt.cols = TRUE, details = FALSE, plot = FALSE)
  expect_equal(round(Test5$`Cmax f2 Factor`$f2,0), 39)
})

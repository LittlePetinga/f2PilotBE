library(testthat)
library(f2PilotBE)


# Test f2.AUC function
test_that("Test for f2.AUC", {

  # Each test validates the computation of f2 factor using different mean
  # concentration-time profiles of Test and Reference products.
  # The expected values are based on predefined expectations or empirical data.

  # Test 1: Verify AUC f2 computation for Conc1
  # Using stacked treatment information
  # Where Reference presents an higher bioavailability
  # But within predefined limits
  Conc1 <- data.frame(
    time = rep(c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
                 2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24), 2),
    trt = c(rep('R', 20), rep('T', 20)),
    conc = c(0.00, 232.44, 389.14, 502.87, 569.75,
             647.44, 688.6, 649.93, 651.57, 661.07,
             605.7, 631.53, 610.52, 570.68, 577.12,
             544.44, 397.27, 308.15, 163.6, 27.49,
             0.00, 202.92, 344, 427.5, 520.84, 559.58,
             553.97, 576.77, 569.9, 582.17, 538.93,
             533.62, 516.38, 481.66, 483.78, 490.83,
             357.08, 261.66, 147.51, 24.06)
  )
  Test1 <- f2.AUC(Conc1, Time = 'time', Conc = "conc",
                  Trt = "trt", Ref = 'R', Test = 'T',
                  Trt.cols = FALSE, details = FALSE, plot = FALSE)
  expect_equal(round(Test1$`AUC f2 Factor`$f2,0), 64)


  # Test 2: Verify AUC f2 computation for Conc2
  # Using stacked treatment information
  # Where Test and Reference presents a similar bioavailability
  Conc2 <- data.frame(
    time = rep(c(0, 0.5, 1, 1.5, 2, 2.5, 2.75, 3, 3.25, 3.5,
                 4, 4.5, 5, 6, 8, 10, 12, 24, 48, 72), 2),
    trt = c(rep('A', 20), rep('B', 20)),
    conc = c(0, 203.21, 350.86, 438.96, 500.61,
             566.3, 574.94, 590.49, 587.42, 567.83,
             547.8, 543.44, 550.92, 533.57, 489.12,
             469.19, 360.07, 266.79, 149.82, 24.02,
             0, 201.93, 331.8, 432.05, 517.49,
             573.08, 570.38, 578.43, 585.06, 587.71,
             588.71, 518.63, 556.32, 507.68, 501.3,
             486.18, 368.36, 273.6, 146.13, 24.08)
  )
  Test2 <- f2.AUC(Conc2, Time = 'time', Conc = "conc",
                  Trt = "trt", Ref = 'A', Test = 'B',
                  Trt.cols = FALSE, details = FALSE, plot = FALSE)
  expect_equal(round(Test2$`AUC f2 Factor`$f2,0), 99)


  # Test 3: Verify AUC f2 computation for Conc3
  # Using stacked treatment information
  # Where Reference presents an higher bioavailability (and not bioequivalent)
  Conc3 <- data.frame(
    Time = rep(c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
                 2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24), 2),
    Trt = c(rep('Ref', 20), rep('Test', 20)),
    Conc = c(0.00, 191.75, 335.47, 414.53, 481.18,
             531.02, 557.02, 555.47, 574.87, 565.03,
             527.63, 494.23, 503.73, 496.94, 469.1,
             445.34, 345.53, 269.26, 140.03, 23.08,
             0.00, 68.87, 126.96, 177.79, 216.37,
             286.97, 318.61, 332.06, 373.68, 387.45,
             388.55, 368.73, 396.64, 406.45, 399.83,
             424.02, 383.47, 331.64, 197.54, 34.84)
  )
  Test3 <- f2.AUC(Conc3, Time = 'Time', Conc = "Conc",
                  Trt = "Trt", Ref = 'Ref', Test = 'Test',
                  Trt.cols = FALSE, details = FALSE, plot = FALSE)
  expect_equal(round(Test3$`AUC f2 Factor`$f2,0), 51)


  # Test 4: Verify AUC f2 computation for Conc4
  # Using stacked treatment information
  # Where Reference presents a lower bioavailability (and not bioequivalent)
  Conc4 <- data.frame(
    Time = rep(c(0, 0.5, 1, 1.5, 2, 2.5, 2.75, 3, 3.25, 3.5,
                 4, 4.5, 5, 6, 8, 10, 12, 24, 48, 72), 2),
    Trt = c(rep('R', 20), rep('T', 20)),
    Conc = c(0.00, 63.15, 116.88, 166.24, 213.65,
             273.23, 287.26, 327.76, 341.85, 361.03,
             392.93, 380.71, 383.28, 390.89, 397.26,
             387.89, 369.19, 315.05, 189.83, 34.39,
             0.00, 191.98, 338.61, 414.11, 470.47,
             551.51, 569.75, 567.95, 585.49, 539.7,
             524.05, 528.5, 523.38, 484.69, 485.84,
             459.27, 352.45, 258.51, 143.19, 23.96)
  )
  Test4 <- f2.AUC(Conc4, Time = 'Time', Conc = "Conc",
                  Trt = "Trt", Ref = 'R', Test = 'T',
                  Trt.cols = FALSE, details = FALSE, plot = FALSE)
  expect_equal(round(Test4$`AUC f2 Factor`$f2,0), 64)


  # Test 5: Verify AUC f2 computation for Conc5
  # Using stacked treatment information
  # Where Reference presents an higher bioavailability (and not bioequivalent)
  Conc5 <- data.frame(
    Time = c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
             2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24),
    Reference = c(0, 221.23, 377.19, 494.73, 555.74,
                  623.86, 615.45, 663.38, 660.29, 621.71,
                  650.33, 622.28, 626.72, 574.94, 610.51,
                  554.02, 409.14, 299.76, 162.85, 27.01),
    Test = c(0, 149.24, 253.05, 354.49, 412.49,
             530.07, 539.68, 566.3, 573.54, 598.33,
             612.63, 567.48, 561.1, 564.47, 541.5,
             536.92, 440.32, 338.78, 185.03, 31.13)
  )
  Test5 <- f2.AUC(Conc5, Time = 'Time',Ref = 'Reference', Test = 'Test',
                  Trt.cols = TRUE, details = FALSE, plot = FALSE)
  expect_equal(round(Test5$`AUC f2 Factor`$f2,0), 71)
})

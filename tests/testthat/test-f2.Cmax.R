test_that("Test for f2.Cmax", {
  dta <- data.frame(
    Time = c(0, 0.25, 0.5, 0.75, 1, 1.5, 1.75, 2, 2.25, 2.5,
             2.75, 3, 3.25, 3.5, 3.75, 4, 6, 8, 12, 24),
    Reference = c(0.00, 221.23, 377.19, 494.73, 555.74, 623.86, 615.45, 663.38, 660.29, 621.71,
                  650.33, 622.28, 626.72, 574.94, 610.51, 554.02, 409.14, 299.76, 162.85, 27.01),
    Test = c(0.00, 149.24, 253.05, 354.49, 412.49, 530.07, 539.68, 566.30, 573.54, 598.33,
             612.63, 567.48, 561.10, 564.47, 541.50, 536.92, 440.32, 338.78, 185.03, 31.13)
  )
  Test_Res <- f2.Cmax(dta, Time = 'Time', Ref = 'Reference', Test = 'Test',
                      Trt.cols = TRUE, details = TRUE, plot = FALSE)
  Test_Res <- round(Test_Res$`Cmax f2 Factor`$f2,4)

  expect_equal(Test_Res, 38.9735)
})
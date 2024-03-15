context("f2")

test_that("Test for f2", {
  dta <- data.frame(
    Time      = c(0, 15, 30, 45, 60, 75, 90),
    Reference = c(0, 40, 67, 80, 87, 89, 91),
    Test1     = c(0, 28, 51, 71, 88, 92, 94),
    Test2     = c(0, 36, 69, 84, 89, 93, 95),
    Test3     = c(0, 43, 78, 86, 93, 94, 96),
    Test4     = c(0, 78, 89, 91, 93, 95, 98)
  )
  Test1 <- f2(dta, Ref = 'Reference', Test = 'Test1')
  Test2 <- f2(dta, Ref = 'Reference', Test = 'Test2')
  Test3 <- f2(dta, Ref = 'Reference', Test = 'Test3')
  Test4 <- f2(dta, Ref = 'Reference', Test = 'Test4')

  expect_equal(round(Test1$f2,0), 52)
  expect_equal(round(Test2$f2,0), 72)
  expect_equal(round(Test3$f2,0), 59)
  expect_equal(round(Test4$f2,0), 36)
})

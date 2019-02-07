context("test-stat")

test_that("testing the proper renaming function", {
  s <- data.frame(q_120_time_between_v4=c("unknown","-9","00 – 10 Minutes"))
  expect_equal(
    s %>% correctLabels1() %>% {.$q_120_time_between_v4[3]},
   factor("10", levels = c("10","30",
                           "60","120","4h","more","4h-5"))
  )
  expect_true(s %>% correctLabels1() %>% {.$q_120_time_between_v4[1]} %>% is.na)
})

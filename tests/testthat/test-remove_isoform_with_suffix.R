# test remove_isoform_suffix
test_that('test_remove_isoform_suffix', {
  dat <- data.frame(a=c(1, 2, 3, 4, 5), b=c(1, 2, 3, 4, 5), row.names = c('P1001', 'P1001-1', 'P1002', 'P1002-2', 'P1003'))
  mat <- as.matrix(dat)
  expected_1 <- as.matrix(data.frame(a=c(1, 3, 5), b=c(1, 3, 5), row.names = c('P1001', 'P1002', 'P1003')))
  actual_1 <- remove_isoform_with_suffix(mat, method='remove')
  diff <-  sum(expected_1 - actual_1)
  expect_equal(diff, 0)
  
  
  
  expected_2 <- as.matrix(data.frame(a=c(3, 7, 5), b=c(3, 7, 5), row.names = c('P1001', 'P1002', 'P1003')))
  actual_2 <- remove_isoform_with_suffix(mat, method='combine')
  # browser()
  diff_2 <- sum(expected_2 - actual_2)
  
  expect_equal(diff_2, 0)
})

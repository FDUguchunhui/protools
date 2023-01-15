testthat::test_that('test_CPM', {
  count <-  matrix(c(10, 12, 30,
                      20, 25, 60,
                      5, 8, 15,
                      0, 0, 1),
                    nrow = 4,
                    byrow = T,
                    dimnames = list(c('A', 'B', 'C', 'D'), c('Rep1', 'Rep2', 'Rep3')))
  expected <- matrix(c(3.3333, 2.9630, 3.3259,
                       3.3333, 3.0864, 3.3259,
                       3.3333, 3.9506, 3.3259,
                       0, 0, 0.02), 
                     byrow = T,
                     nrow = 4,
                     dimnames = list(c('A', 'B', 'C', 'D'), c('Rep1', 'Rep2', 'Rep3')))
  gene_length <- data.frame(length=c(10, 1, 4, 2), row.names = c('D', 'C', 'B', 'A'))
  actual <- TPM(count, gene_length=gene_length, per_count = 10)
  diff <-  sum(actual - expected)
  testthat::expect_equal(diff, 0, tolerance = 0.01)
})



# test mapping function
testthat::test_that('test_mapping', {
  dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
  keys <- c('a', 'a', 'c', 'b')
  expected <- c('A', 'A', 'C', 'B')
  actual <- mapping(keys, dictionary)
  testthat::expect_equal(expected, actual)
  
  # when dictionary has non-unique mapping for a key
  dictionary_2 <- tibble::tibble(key=c('a', 'a', 'b', 'c'), value=c('A', 'B',  'B', 'C'))
  testthat::expect_error(mapping(keys, dictionary_2))
  
})


# test column_mapping function
testthat::test_that('test_column_mapping', {
  dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
  A <- tibble::tibble(letter=c('a', 'a', 'c', 'b'), noun=c('apple', 'apple', 'cab', 'bean'))
  actual <- column_mapping(A, 'letter', dictionary)
  expected <- tibble::tibble(letter=c('A', 'A', 'C', 'B'), noun=c('apple', 'apple', 'cab', 'bean'))
  testthat::expect_equal(actual, expected)
})



# test remove_isoform_suffix
testthat::test_that('test_remove_isoform_suffix', {
  dat <- data.frame(a=c(1, 2, 3, 4, 5), b=c(1, 2, 3, 4, 5), row.names = c('P1001', 'P1001-1', 'P1002', 'P1002-2', 'P1003'))
  mat <- as.matrix(dat)
  expected_1 <- as.matrix(data.frame(a=c(1, 3, 5), b=c(1, 3, 5), row.names = c('P1001', 'P1002', 'P1003')))
  actual_1 <- remove_isoform_suffix(mat, method='remove')
  diff <-  sum(expected_1 - actual_1)
  testthat::expect_equal(diff, 0)
  
  
  
  expected_2 <- as.matrix(data.frame(a=c(3, 7, 5), b=c(3, 7, 5), row.names = c('P1001', 'P1002', 'P1003')))
  actual_2 <- remove_isoform_suffix(mat, method='combine')
  browser()
  diff_2 <- sum(expected_2 - actual_2)

  testthat::expect_equal(diff_2, 0)
})



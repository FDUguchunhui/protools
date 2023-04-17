# test column_mapping function
test_that('test_column_mapping', {
  dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
  A <- tibble::tibble(letter=c('a', 'a', 'c', 'b'), noun=c('apple', 'apple', 'cab', 'bean'))
  actual <- column_mapping(A, 'letter', dictionary)
  expected <- tibble::tibble(letter=c('A', 'A', 'C', 'B'), noun=c('apple', 'apple', 'cab', 'bean'))
  expect_equal(actual, expected)
})




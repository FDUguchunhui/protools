test_that("mapping 'keys' to 'values' based on dictionary and return with the 
          same order as in 'keys'", {
  # browser()
  dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
  keys <- c('a', 'a', 'c', 'b')
  actual <- mapping(keys, dictionary)
  expected <- c('A', 'A', 'C', 'B')
  expect_equal(actual, expected)
})


test_that("when an empty 'keys' vector is used", {
            dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
            keys <- character()
            actual <- mapping(keys, dictionary)
            expected <- character()
            expect_equal(actual, expected)
          })


test_that("when duplicated keys exists in the dictionary", {
  dictionary <- tibble::tibble(key=c('a', 'a', 'c'), value=c('A', 'B', 'C'))
  keys <- c('a', 'a', 'c', 'b')
  expect_error(mapping(keys, dictionary))
})


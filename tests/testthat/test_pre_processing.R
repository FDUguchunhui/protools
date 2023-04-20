# test_that('test_CPM', {
#   count <-  matrix(c(10, 12, 30,
#                       20, 25, 60,
#                       5, 8, 15,
#                       0, 0, 1),
#                     nrow = 4,
#                     byrow = T,
#                     dimnames = list(c('A', 'B', 'C', 'D'), c('Rep1', 'Rep2', 'Rep3')))
#   expected <- matrix(c(3.3333, 2.9630, 3.3259,
#                        3.3333, 3.0864, 3.3259,
#                        3.3333, 3.9506, 3.3259,
#                        0, 0, 0.02), 
#                      byrow = T,
#                      nrow = 4,
#                      dimnames = list(c('A', 'B', 'C', 'D'), c('Rep1', 'Rep2', 'Rep3')))
#   gene_length <- data.frame(length=c(10, 1, 4, 2), row.names = c('D', 'C', 'B', 'A'))
#   actual <- TPM(count, gene_length=gene_length, per_count = 10)
#   diff <-  sum(actual - expected)
#   testthat::expect_equal(diff, 0, tolerance = 0.01)
# })
# 



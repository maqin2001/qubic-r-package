library(QUBIC)
context("Sample")

expect_equal(10, 10)
expect_equal(10, 10 + 1e-7)

test <- matrix(rnorm(5000),100,50)
test[11:20,11:20] <- rnorm(100,3,0.3)
res <- biclust::biclust(test, method=BCQU())
expect_true(isS4(res))
expect_true(is.object(res))
res

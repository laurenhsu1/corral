context("corral")

library(corral)
library(ade4)

test_that('same eigens as dudi.coa',{
  mat <- matrix(sample(0:10, 500, replace=TRUE), ncol=25)
  a <- corral(mat,ncomp = 2)
  b <- dudi.coa(mat, scannf = FALSE, nf = 2)
  expect_equal(all.equal(b$eig[1:2],a$d[1:2]^2),TRUE)
})

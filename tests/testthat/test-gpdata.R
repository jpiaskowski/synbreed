context('gpData object type and associated functions')

set.seed(123)
# 9 plants with 2 traits
n <- 9 # only for n > 6
pheno <- data.frame(Yield = rnorm(n, 200, 5), Height = rnorm(n, 100, 1))
rownames(pheno) <- letters[1:n]

# marker matrix
geno <- matrix(sample(c("AA", "AB", "BB", NA),
  size = n * 12, replace = TRUE,
  prob = c(0.6, 0.2, 0.1, 0.1)), nrow = n)

rownames(geno) <- letters[n:1]
colnames(geno) <- paste("M", 1:12, sep = "")

# genetic map
# one SNP is not mapped (M5) and will therefore be removed
map <- data.frame(chr = rep(1:3, each = 4), pos = rep(1:12))
map <- map[-5, ]
rownames(map) <- paste("M", c(1:4, 6:12), sep = "")

# simulate pedigree
ped <- simul.pedigree(3, c(3, 3, n - 6))

# combine in one object
gp <- create.gpData(pheno, geno, map, ped, cores = 2)
summary(gp)

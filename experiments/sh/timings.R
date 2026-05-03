set.seed(42)
library(balnet)
library(microbenchmark)
library(xtable)

headers = names(read.csv("treated2008_conifer_X.csv.gz", nrows = 1))
X = matrix(
  scan("treated2008_conifer_X.csv.gz", sep = ",", skip = 1),
  ncol = length(headers),
  dimnames = list(NULL, headers),
  byrow = TRUE
)
W = scan("treated2008_conifer_W.csv.gz", skip = 1)

get_data = function(size = c("100k", "200k", "400k")) {
  N = 100000
  P = 500
  if (size == "100k") {
    W.out = W[1:N]
    X.out = X[1:N, 1:P]
  } else if (size == "200k") {
    n = 2 * N
    row = sample.int(N, n, replace = TRUE)
    col = rep(1:P, 2)
    W.out = W[row]
    X.out = X[row, col]
  } else if (size == "400k") {
    n = 4 * N
    row = sample.int(N, n, replace = TRUE)
    col = rep(1:P, 4)
    W.out = W[row]
    X.out = X[row, col]
  }
  colnames(X.out) <- make.names(1:ncol(X.out))

  list(W = W.out, X = X.out)
}

times = 2
num.threads = 4
n = setNames(c(1e5, 2*1e5, 4*1e5), c("100k", "200k", "400k"))
p = setNames(c(500, 2*500, 4*500), c("100k", "200k", "400k"))
grid = expand.grid(
  size = c("100k", "200k", "400k"),
  max.imbalance = c(0.05, 0.01),
  stringsAsFactors = FALSE
)
out = list()
for (i in 1:nrow(grid)) {
  size = grid$size[i]
  max.imbalance = grid$max.imbalance[i]
  data = get_data(size = size)
  bench = microbenchmark(
    balnet(data$X, data$W, target = "ATT", max.imbalance = max.imbalance, num.threads = num.threads),
    times = times,
    unit = "seconds"
  )
  out = c(out, list(data.frame(n = n[[size]], p = p[[size]], max.imbalance = max.imbalance, time = summary(bench)$mean)))
}
out.df = do.call(rbind, out)
write.csv(out.df, "timing.csv", row.names = FALSE)

# Raw numbers in seconds behind Table 1
print(xtable(out.df, digits = 1), include.rownames = FALSE)

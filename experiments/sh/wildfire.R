library(balnet)

colnames = names(read.csv("treated2008_conifer_X.csv.gz", nrows = 1))
X = matrix(
  scan("treated2008_conifer_X.csv.gz", sep = ",", skip = 1),
  ncol = length(colnames),
  dimnames = list(NULL, colnames),
  byrow = TRUE
)
W = scan("treated2008_conifer_W.csv.gz", skip = 1)

# Maximum unweighted imbalance is simply lambda_max
print(balnet(X, W, target = "ATT", nlambda = 1))

# Fit a balnet targeting the ATT.
fit <- balnet(X, W, target = "ATT", max.imbalance = 0.05, verbose = TRUE, num.threads = 4)
print(fit)

# Plot path
pdf("fig_path.pdf", pointsize = 18)
plot(fit)
dev.off()

# Plot SMD at lowest lambda
pdf("fig_smd_top.pdf", pointsize = 18)
par(mar = c(5, 5.5, 4, 2))
plot(fit, lambda = 0, max = 10)
dev.off()

# Plot all covariates, grouped
covariate.groups = list(
  `Min air temperature` = grep("minat_", colnames(X)),
  `Max air temperature` = grep("maxat_", colnames(X)),
  `Precipitation` = grep("prcp_", colnames(X)),
  `Snow water equivalent` = grep("swe_", colnames(X)),
  `Water vapor pressure` = grep("wvp_", colnames(X)),
  `Fire frequency` = grep("fire_\\d", colnames(X)),
  `Avg fire brightness` = grep("avg_BRIGHTNESS_", colnames(X)),
  `Max fire radiative power` = grep("max_FRP_", colnames(X)),
  `Disturbance: fire` = grep("fire_disturb_", colnames(X)),
  `Disturbance: timber harvest` = grep("timber_", colnames(X)),
  `Disturbance: drought` = grep("drought_", colnames(X)),
  `Disturbance: greening` = grep("greening_", colnames(X)),
  `Disturbance: browning` = grep("browning_", colnames(X)),
  `Vegetation: tree cover` = grep("tree_cover_", colnames(X)),
  `Elevation` = grep("elev", colnames(X))
)
pdf("fig_smd_grp.pdf", pointsize = 18)
par(mar = c(5, 8.5, 4, 2))
plot(fit, lambda = 0, groups = covariate.groups)
dev.off()

# Inspect ATT weights at min. lambda
weights = balweights(fit, lambda = 0)

png("fig_weights.png", width = 2100, height = 2100, res = 300, pointsize = 18)
hist(weights$control[W==0], breaks = 150, prob = TRUE,
     col = rgb(0,0,0,0.4),
     main = "",
     xlab = "")
rug(weights$control[W==0], col = rgb(0,0,0,0.4))
dev.off()

# cv.balnet example
set.seed(42)
cv.fit <-  cv.balnet(X, W, target = "ATT", nfolds = 5, max.imbalance = 0.05, verbose = TRUE, num.threads = 4)
print(cv.fit)

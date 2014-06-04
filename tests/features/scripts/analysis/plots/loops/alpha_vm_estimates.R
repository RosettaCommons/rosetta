library(RSQLite)

db_path <- "/scratch/weitzner/loop_features_test/antibodies_f33ffc7_130904/features_antibodies_f33ffc7_130904.db3"

sele <-paste("SELECT alpha101 FROM loop_anchor_transforms")


sqlite    <- dbDriver("SQLite")
db <- dbConnect(sqlite, db_path)

f <- dbGetQuery(db, sele)

dbDisconnect(db)

library(movMF)
#movMF(data.matrix(f$alpha101), 2)

#library(mixtools)
library(ggthemes)
# f$alpha101_scaled <- ifelse(f$alpha101 < -100, f$alpha101 + 360, f$alpha101)
# 
# expectation_maximization.model <- normalmixEM(x=f$alpha101_scaled)
# parameters <- expectation_maximization.model[c("lambda", "mu", "sigma")]
# 
# f$curve.1 <- parameters$lambda[1] * dnorm(f$alpha101_scaled, parameters$mu[1],
#                                           parameters$sigma[1])
# 
# f$curve.2 <- parameters$lambda[2] * dnorm(f$alpha101_scaled, parameters$mu[2],
#                                           parameters$sigma[2])
# 
# f$combined <- f$curve.1 + f$curve.2

pts_on_unit_circle <- cbind(cos(f$alpha101 * pi / 180), 
                            sin(f$alpha101 * pi / 180))

d <- movMF(pts_on_unit_circle, 2)

# detect quadrant of theta for each component of mixture
norm_theta <- skmeans:::row_normalize(d$theta)

# compute mean angles of distributions
mu <- atan2(norm_theta[,2], norm_theta[,1]) * 180 / pi

# compute the standard deviations
kappa <- (d$theta / norm_theta)[,1]


# Because the data are laid out on a unit circle and we will be plotting
# on a domain of 360, divide the densities by 360.
f$curve.1 <- d$alpha[1] * dmovMF(pts_on_unit_circle, d$theta[1,]) / 360
f$curve.2 <- d$alpha[2] * dmovMF(pts_on_unit_circle, d$theta[2,]) / 360
f$combined <- f$curve.1 + f$curve.2

hline.data <- data.frame(z = seq(0.005, 0.025, by=0.005))
p <- ggplot(f, aes(x=alpha101))
p + geom_histogram(aes(y=..density..), fill="lightgrey", colour="white") + 
#   geom_line(aes(y=combined), size=1, colour="black") +
  geom_hline(aes(yintercept = z), hline.data, colour="white") +
  geom_line(aes(y=curve.1), size=1, colour="grey50") +
  geom_line(aes(y=curve.2), size=1, colour="darkgrey") + 
  theme_tufte() + #geom_rangeframe() +
  scale_x_continuous(expression(paste(alpha[101], " (degrees)")), 
                     limits=c(-180, 180), breaks=seq(-180,180, by=45)) +
  scale_y_continuous("Density", expand=c(0, 0))

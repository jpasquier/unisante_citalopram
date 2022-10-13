library(survival)

# number of observations
N <- 10^3

# Mean and SD of survival
m <- 50
s <- 10

# Not 'left truncated' data
T0 <- rnorm(N, m, s)
S0 <- Surv(rep(0, length(T0)), T0, rep(1, length(T0)))

# Left truncated data
T1 <- rnorm(2 * N, m, s)
T1_start <- runif(length(T1), min(T1), max(T1))
b  <- T1 >= T1_start
S1 <- Surv(T1_start[b], T1[b], rep(1, sum(b)))

# Comparison
dta <- data.frame(
  S = c(S0, S1),
  X = c(rep(0, length(S0)), rep(1, length(S1)))
)
fit <- survfit(S ~ X, dta)
x11();plot(fit)

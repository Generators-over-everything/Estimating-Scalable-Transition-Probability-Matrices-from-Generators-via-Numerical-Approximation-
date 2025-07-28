###############################################################
# 1) load packages (assumes you've already run install.packages)
###############################################################
library(markovMSM)
library(mstate)
library(survival)

###############################################################
# 2) load ebmt4 and build tmat exactly as before
###############################################################
data("ebmt4")
db_wide <- ebmt4

positions <- list(
  c(2, 3, 5, 6),
  c(4, 5, 6),
  c(4, 5, 6),
  c(5, 6),
  c(6),
  c()
)
state.names <- c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
tmat <- transMatMSM(positions, state.names)

timesNames <- c(NA, "rec", "ae", "recae", "rel", "srv")
status    <- c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s")

###############################################################
# 3) convert to “long” format
###############################################################
db_long <- prepMSM(
  data       = db_wide,
  trans      = tmat,
  timesNames = timesNames,
  status     = status
)

# — remove any rows where Tstop <= Tstart so that Surv(...) is valid:
db_long <- subset(db_long, Tstop > Tstart)

# ensure “trans” is a factor (msfit/coxph prefer that)
db_long$trans <- as.factor(db_long$trans)

###############################################################
# 4) fit the Cox stratified by transition
###############################################################
cox_fit <- coxph(
  Surv(Tstart, Tstop, status) ~ strata(trans),
  data = db_long
)

###############################################################
# 5) pass “trans = tmat” into msfit() so it knows the structure
###############################################################
cumhaz_list <- msfit(cox_fit, trans = tmat)

# cumhaz_list now contains a cumulative‐hazard object for each allowed transition

###############################################################
# 6) build a constant‐rate 6×6 Q matrix from those cumhazards
###############################################################
library(dplyr)

# extract cumulative hazards via basehaz()
bh_all <- basehaz(cox_fit, centered = FALSE)
bh_split <- split(bh_all, bh_all$strata)

Qmat <- matrix(0, nrow = 6, ncol = 6)
rownames(Qmat) <- colnames(Qmat) <- state.names

for (j in seq_along(bh_split)) {
  df_j <- bh_split[[j]]
  last_row <- df_j[nrow(df_j), ]
  t_max <- last_row$time
  H_j   <- last_row$hazard
  λ_j   <- H_j / t_max           # constant‐rate approximation

  idx <- which(tmat == j, arr.ind = TRUE)
  i   <- idx[1]
  k   <- idx[2]
  Qmat[i, k] <- λ_j
}

for (i in 1:6) {
  Qmat[i, i] <- -sum(Qmat[i, -i])
}

print(Qmat)

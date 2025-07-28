# -------------------------------
# 1. Package Installation / Loading
# -------------------------------

library(expm)
library(Matrix)
library(xlsx)

# -------------------------------
# 2. Sonnenberg Method
# -------------------------------

Sonnenberg <- function(A, p) {
  n <- dim(A)[1]
  ii <- 1:n

  for (i in ii) {
    for (j in ii) {
      if (j != i) {
        A[i, j] <- 1 - (1 - A[i, j])^p
      }
    }
    A[i, i] <- 1 - sum(A[i, (1:n)]) + A[i, i]
  }

  return(A)
}

# -------------------------------
# 3. Schur–Padé Method (SchurPade)
# -------------------------------

SchurPade <- function(A, p) {
  maxsqrt <- 36
  n      <- dim(A)[1]
  nsq    <- 0
  m      <- 0

  # Compute (complex) Schur decomposition A = Q T Q*
  AS    <- Schur(A)
  Q     <- AS$Q
  Tmat  <- AS$T
  diagT <- diag(Tmat)

  # Check eigenvalues of A
  if (nnzero(diagT) != n) {
    warning("Matrix power may not exist for singular matrix")
  }
  if (any(Im(diagT) == 0 & Re(diagT) <= 0)) {
    warning("Principal matrix power is not defined for A with nonpositive eigenvalues. A non-principal power is returned")
  }

  # If T is already diagonal, just do diag powers
  if (identical(Tmat, diag(diagT, n))) {
    X <- Q %*% diag(diagT^p, n) %*% t(Conj(Q))
    return(X)
  }

  #############################
  # 3.1. Helper: powerm2by2
  #############################
  powerm2by2 <- function(A2, p2) {
    a1  <- as.complex(A2[1, 1])
    a2  <- as.complex(A2[2, 2])
    a1p <- a1^p2
    a2p <- a2^p2
    loga1 <- log(a1)
    loga2 <- log(a2)

    X <- diag(c(a1p, a2p))

    if (a1 == a2) {
      X[1, 2] <- p2 * A2[1, 2] * a1^(p2 - 1)
      return(X)
    } else if (abs(a1) < 0.5 * abs(a2) || abs(a2) < 0.5 * abs(a1)) {
      X[1, 2] <- A2[1, 2] * (a2p - a1p) / (a2 - a1)
      return(X)
    } else {
      #########################
      # 3.1.1. Helper: unwinding
      #########################
      unwinding <- function(z, k = 0) {
        if (k == 0) {
          u <- ceiling((Im(z) - pi) / (2 * pi))
          return(u)
        } else {
          return(0)
        }
      }

      w  <- atanh((a2 - a1) / (a2 + a1)) + 1i * pi * unwinding(loga2 - loga1, 0)
      dd <- 2 * exp(p2 * (loga1 + loga2) / 2) * sinh(p2 * w) / (a2 - a1)
      X[1, 2] <- A2[1, 2] * dd
      return(X)
    }
  }

  #############################
  # 3.2. Helper: sqrtm_tri
  #############################
  sqrtm_tri <- function(Ttri) {
    n <- dim(Ttri)[1]
    # Convert to complex form
    for (i in 1:n) {
      for (j in 1:n) {
        Ttri[i, j] <- as.complex(Ttri[i, j])
      }
    }
    R <- diag(0 + 0i, n)

    for (j in 1:n) {
      R[j, j] <- sqrt(Ttri[j, j])
      if (j > 1) {
        for (i in seq(j - 1, 1, -1)) {
          numerator   <- Ttri[i, j] - R[i, (i + 1):(j - 1)] %*% R[(i + 1):(j - 1), j]
          denominator <- R[i, i] + R[j, j]
          R[i, j]     <- numerator / denominator
        }
      }
    }
    return(R)
  }

  #############################
  # 3.3. Helper: coeff
  #############################
  coeff <- function(pcoeff, icoeff) {
    if (icoeff == 1) {
      return(-pcoeff)
    } else {
      g <- icoeff / 2
      if (g == round(g)) {
        cval <- (-g + pcoeff) / (2 * (2 * g - 1))
        return(cval)
      } else {
        g <- floor(g)
        cval <- (-g - pcoeff) / (2 * (2 * g + 1))
        return(cval)
      }
    }
  }

  #############################
  # 3.4. Helper: powerm_cf
  #############################
  powerm_cf <- function(Y, pcf, mcf) {
    # Evaluate [mcf/mcf] Pade approximant of (I - Y)^p via continued fraction
    k <- 2 * mcf
    n <- dim(Y)[1]

    S <- coeff(pcf, k) * Y
    d <- (k - 1):1
    for (i in d) {
      z <- coeff(pcf, i)
      S <- z * solve((diag(1, n) + S)) %*% Y
    }
    S <- diag(1, n) + S
    return(S)
  }

  #############################
  # 3.5. Helper: powerm_triang (fixed)
  #############################
  powerm_triang <- function(Ttri, ptri, maxsqrt) {
    n <- dim(Ttri)[1]
    nsq <- 0

    if (n == 1) {
      return(Ttri^ptri)
    }
    if (n == 2) {
      return(powerm2by2(Ttri, ptri))
    }

    T_old <- Ttri
    xvals <- c(
      1.512666672122460e-005, 2.236550782529778e-003, 1.882832775783885e-002,
      6.036100693089764e-002, 1.239372725584911e-001, 1.998030690604271e-001,
      2.787629930862099e-001, 3.547373395551596e-001, 4.245558801949280e-001,
      4.870185637611313e-001, 5.420549053918690e-001, 5.901583155235642e-001,
      6.320530128774397e-001, 6.685149002867240e-001, 7.002836650662422e-001,
      7.280253837034645e-001, 9.152924199170567e-001, 9.764341682154458e-001
    )

    while (TRUE) {
      Mmat     <- as.matrix(Ttri - diag(1, n))
      normdiff <- norm(Mmat, "1")

      if (normdiff <= xvals[7]) {
        j1 <- which(xvals[3:7] >= normdiff)[1] + 2
        j2 <- which(xvals[3:7] >= normdiff / 2)[1] + 2

        # FIXED: no assignment inside "if"
        if (j1 - j2 <= 1) {
          m <- j1
          break
        }
      }

      if (nsq == maxsqrt) {
        m <- 16
        break
      }

      Ttri <- sqrtm_tri(Ttri)
      nsq  <- nsq + 1
    }

    X <- powerm_cf(diag(1, n) - Ttri, ptri, m)
    # Squaring phase
    for (s in 1:nsq) {
      X <- X %*% X
    }
    # Replace 2×2 diagonal blocks
    for (i in 1:(n - 1)) {
      Tii <- T_old[i:(i + 1), i:(i + 1)]
      Si  <- powerm2by2(Tii, ptri / (2^(nsq)))
      X[i:(i + 1), i:(i + 1)] <- Si
    }
    return(X)
  }

  #############################
  # 3.6. Final evaluation in SchurPade
  #############################
  Xtri <- powerm_triang(Tmat, p, maxsqrt)
  X    <- Q %*% Xtri %*% t(Conj(Q))
  return(X)
}

# -------------------------------
# 4. Eigenvalue Method (with tryCatch)
# -------------------------------

Eigenvalue <- function(A, p) {
  n <- dim(A)[1]
  # Attempt to compute logm; if it errors, return an NA n×n matrix
  G <- tryCatch({
    logm(A, method = "Eigen")
  }, error = function(e) {
    warning("Eigenvalue::logm failed (non‐diagonalizable or eigenvalue issue); returning NA matrix")
    return(matrix(NA_real_, n, n))
  })

  # If G is all NA, propagate NA
  if (all(is.na(G))) {
    return(matrix(NA_real_, n, n))
  }

  # Otherwise proceed to exponentiate
  E <- expm(G * p)
  return(E)
}

# -------------------------------
# 5. Relative Entropy Regularization (RegRE)
# -------------------------------

RegRE <- function(P) {
  d <- dim(P)[1]
  P <- Re(P)
  z <- which(P < 0)
  P[z] <- 0

  for (i in 1:d) {
    row_sum   <- sum(P[i, 1:d])
    P[i, 1:d] <- P[i, 1:d] / row_sum
  }
  return(P)
}

# -------------------------------
# 6. Extreme Relative Entropy Regularization (RegERE)
# -------------------------------

RegERE <- function(P) {
  d <- dim(P)[1]

  for (i in 1:d) {
    # Compare only the real part
    if ((sum(Re(P[i, 1:d])) - 1) < 1e-16 && min(Re(P[i, 1:d])) >= 0) {
      break
    } else {
      z <- which(Re(P[i, 1:d]) < 0)
      P[i, z] <- 0

      maxP <- max(Re(P[i, 1:d]))
      w    <- which(Re(P[i, 1:d]) != 0)
      minP <- min(Re(P[i, w]))

      if ((maxP + minP > 1) || (maxP == minP)) {
        warning("max + min > 1 or min = max; this algorithm cannot be applied")
        break
      }

      delta <- 1 - maxP - minP
      nsum  <- sum(Re(P[i, 1:d])) - maxP - minP

      for (j in 1:d) {
        if (Re(P[i, j]) != maxP && Re(P[i, j]) != minP) {
          P[i, j] <- (P[i, j] / nsum) * delta
        }
      }
    }
  }
  return(P)
}

# -------------------------------
# 7. Q-Matrix Regularization (RegQ)
# -------------------------------

RegQ <- function(A, p) {
  R   <- logm(A, method = "Eigen")
  d   <- dim(R)[1]
  g   <- matrix(0, d, 1)
  b   <- matrix(0, d, 1)

  for (i in 1:d) {
    for (j in 1:d) {
      if (j != i) {
        g[i] <- g[i] + max(Re(R[i, j]), 0)
        b[i] <- b[i] + max(-Re(R[i, j]), 0)
      }
    }
  }
  for (i in 1:d) {
    g[i] <- abs(Re(R[i, i])) + g[i]
  }

  for (i in 1:d) {
    for (j in 1:d) {
      if (j != i && Re(R[i, j]) < 0) {
        R[i, j] <- 0
      }
      if (g[i] > 0) {
        R[i, j] <- Re(R[i, j]) - b[i] * abs(Re(R[i, j])) / g[i]
      }
      # if g[i] == 0, R[i,j] stays as is
    }
  }

  A1 <- expm(R * p)
  return(A1)
}

# -------------------------------
# 8. Complete convertTP Function (Eigenvalue branch wrapped in tryCatch)
#       (Tracks NA reasons)
# -------------------------------
convertTP <- function(A, p) {
  # (Assumes packages are already installed/loaded.)

  # 8.1. Compute M1, M2 first
  M1 <- Sonnenberg(A, p)
  M2 <- SchurPade(A, p)

  # 8.1 (continued): compute M3 but catch NA case
  M3_candidate <- Eigenvalue(A, p)
  if (all(is.na(M3_candidate))) {
    M3 <- matrix(NA_real_, nrow(A), ncol(A))
  } else {
    M3 <- Re(M3_candidate)
  }

  # Ensure we also take only Re() of M2
  M2 <- Re(M2)

  d <- dim(A)[1]

  # Prepare storage for method errors and candidate matrices
  error    <- matrix(NA_real_, 8, 2)
  error[, 1] <- c("M1", "M2", "M21", "M22", "M3", "M31", "M32", "M33")

  # Candidate matrices (initialize as NA matrices to fill in)
  M21 <- matrix(NA_real_, d, d)
  M22 <- matrix(NA_real_, d, d)
  M31 <- matrix(NA_real_, d, d)
  M32 <- matrix(NA_real_, d, d)
  M33 <- matrix(NA_real_, d, d)

  # 8.2. Check stochastic and compute raw errors for M1
  if (all(is.na(M1))) {
    reason1 <- "failed"
    error[1, 2] <- NA
  } else if (abs(sum(Re(M1)) - d) < 1e-16 && min(Re(M1)) >= 0) {
    reason1 <- "valid"
    error[1, 2] <- norm(A - M1^(1/p), "f") / norm(A, "f")
  } else {
    reason1 <- "non-stochastic"
    error[1, 2] <- NA
  }

  # 8.2. Schur-Padé part for M2
  if (all(is.na(M2))) {
    reason2 <- "failed"
    error[2, 2] <- NA
  } else if (abs(sum(Re(M2)) - d) < 1e-16 && min(Re(M2)) >= 0) {
    reason2 <- "valid"
    error[2, 2] <- norm(A - M2^(1/p), "f") / norm(A, "f")
  } else {
    reason2 <- "non-stochastic"
    error[2, 2] <- NA
  }

  # 8.2. Regularizations for M2 if needed
  if (reason2 == "non-stochastic") {
    M21 <- RegRE(M2)
    if (all(is.na(M21))) {
      reason3 <- "failed"
      error[3, 2] <- NA
    } else if (abs(sum(Re(M21)) - d) < 1e-16 && min(Re(M21)) >= 0) {
      reason3 <- "valid"
      error[3, 2] <- norm(A - M21^(1/p), "f") / norm(A, "f")
    } else {
      reason3 <- "non-stochastic"
      error[3, 2] <- NA
    }

    M22 <- RegERE(M2)
    if (all(is.na(M22))) {
      reason4 <- "failed"
      error[4, 2] <- NA
    } else if (abs(sum(Re(M22)) - d) < 1e-16 && min(Re(M22)) >= 0) {
      reason4 <- "valid"
      error[4, 2] <- norm(A - M22^(1/p), "f") / norm(A, "f")
    } else {
      reason4 <- "non-stochastic"
      error[4, 2] <- NA
    }
  } else {
    reason3 <- "not-applicable"
    reason4 <- "not-applicable"
  }

  # 8.3. Eigenvalue branch for M3
  if (all(is.na(M3_candidate))) {
    reason5 <- "failed"
    error[5, 2] <- NA
  } else if (abs(sum(Re(M3)) - d) < 1e-16 && min(Re(M3)) >= 0) {
    reason5 <- "valid"
    error[5, 2] <- norm(A - M3^(1/p), "f") / norm(A, "f")
  } else {
    reason5 <- "non-stochastic"
    error[5, 2] <- NA
  }

  # 8.3. Regularizations for M3 if needed
  if (reason5 == "non-stochastic") {
    M31 <- RegRE(M3)
    if (all(is.na(M31))) {
      reason6 <- "failed"
      error[6, 2] <- NA
    } else if (abs(sum(Re(M31)) - d) < 1e-16 && min(Re(M31)) >= 0) {
      reason6 <- "valid"
      error[6, 2] <- norm(A - M31^(1/p), "f") / norm(A, "f")
    } else {
      reason6 <- "non-stochastic"
      error[6, 2] <- NA
    }

    M32 <- RegERE(M3)
    if (all(is.na(M32))) {
      reason7 <- "failed"
      error[7, 2] <- NA
    } else if (abs(sum(Re(M32)) - d) < 1e-16 && min(Re(M32)) >= 0) {
      reason7 <- "valid"
      error[7, 2] <- norm(A - M32^(1/p), "f") / norm(A, "f")
    } else {
      reason7 <- "non-stochastic"
      error[7, 2] <- NA
    }

    M33 <- RegQ(A, p)
    if (all(is.na(M33))) {
      reason8 <- "failed"
      error[8, 2] <- NA
    } else if (abs(sum(Re(M33)) - d) < 1e-16 && min(Re(M33)) >= 0) {
      reason8 <- "valid"
      error[8, 2] <- norm(A - M33^(1/p), "f") / norm(A, "f")
    } else {
      reason8 <- "non-stochastic"
      error[8, 2] <- NA
    }
  } else {
    reason6 <- "not-applicable"
    reason7 <- "not-applicable"
    reason8 <- "not-applicable"
  }

  # 8.4. Compile status reasons for each row
  status <- c(
    reason1, reason2, reason3, reason4,
    reason5, reason6, reason7, reason8
  )

  # 8.5. Print the error matrix and status
  colnames(error) <- c("Method", "Error")
  rownames(error) <- c("M1", "M2", "M21", "M22", "M3", "M31", "M32", "M33")
  print(error)

  status_df <- data.frame(
    Method = c("M1", "M2", "M21", "M22", "M3", "M31", "M32", "M33"),
    Status = status,
    stringsAsFactors = FALSE
  )
  print(status_df)

  # 8.6. Find best valid method by smallest error
  valid_rows <- which(status == "valid")
  if (length(valid_rows) == 0) {
    stop("No valid method found for this matrix.")
  }
  best_idx      <- valid_rows[ which.min(error[valid_rows, 2]) ]
  chosen_method <- status_df$Method[best_idx]
  cat("Chosen method:", chosen_method, "\n")

  # 8.7. Return the chosen method name and its matrix
  Solution <- get(chosen_method)
  return(list(
    method = chosen_method,
    matrix = Solution
  ))
}

# -------------------------------
# 9. Example Usage
# -------------------------------

P_test <- matrix(
  c(
    0.20, 0.30, 0.10, 0.20, 0.20,
    0.10, 0.40, 0.20, 0.20, 0.10,
    0.25, 0.25, 0.25, 0.15, 0.10,
    0.30, 0.20, 0.20, 0.20, 0.10,
    0.00, 0.00, 0.00, 0.00, 1.00
  ),
  byrow = TRUE,
  nrow = 5
)

P_monthly <- convertTP(P_test, 12)
print("Original P_test:")
print(P_test)
print("Returned p-cycle matrix and method:")
print(P_monthly)

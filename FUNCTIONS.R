### Hamilton - Regime Switching Markov Model

# matf(): construye la matriz de transición F (n x n) del estado expandido
# a partir de la matriz de transición primitiva pm (ns x ns),
# siguiendo exactamente el algoritmo iterativo del GAUSS original.
#
# Args:
#   pm : matriz ns x ns con probabilidades de transición primitivas
#   ns : número de estados primitivos
#   ps : número de rezagos del estado en el vector de estado (ps >= 0)
#
# Return:
#   fm : matriz F de tamaño n x n, con n = ns^(ps+1)

# ====================== PROCS.R (solo funciones) ======================

# ---- matpm: parámetros -> matriz de transición pm (ns x ns) ----
matpm <- function(xth, ns, ipm, izz) {
  tr <- function(x) if (izz == 1L) x^2/(1 + x^2) else x
  pm <- matrix(0, ns, ns)
  
  if (ipm == 1L) {
    stopifnot(ns == 2L, length(xth) >= 2L)
    p11 <- tr(xth[1]); p22 <- tr(xth[2])
    pm[1,1] <- p11;           pm[2,2] <- p22
    pm[2,1] <- 1 - p11;       pm[1,2] <- 1 - p22
    return(pm)
  }
  
  if (ipm == 2L) {
    stopifnot(length(xth) == (ns - 1L) * ns)
    pm[1:(ns - 1L), ] <- matrix(xth, nrow = ns - 1L, ncol = ns, byrow = FALSE)
    if (izz == 1L) {
      pm[ns, ] <- 1
      pm <- pm^2
      pm <- sweep(pm, 2, colSums(pm), "/")  # normaliza columnas
    } else {
      pm[ns, ] <- 1 - colSums(pm)
    }
    return(pm)
  }
  
  if (ipm == 3L) {
    stopifnot(ns == 3L, length(xth) == 3L)
    p11 <- tr(xth[1]); p12 <- tr(xth[2]); p23 <- tr(xth[3])
    pm[1, ] <- c(p11, 1 - p11, 0)
    pm[2, ] <- c(0,    1 - p12, p23)
    pm[3, ] <- c(1 - p11, 0,    1 - p23)  
    return(pm)
  }
  
  stop("ipm no soportado en matpm().")
}

# ---- matf: matriz de transición en el estado expandido ----
matf <- function(pm, ns, ps) {
  stopifnot(is.matrix(pm), nrow(pm) == ns, ncol(pm) == ns, ns >= 2, ps >= 0)
  na <- 1L; nb <- ns; nc <- ns * ns; fm <- pm
  if (ps == 0L) return(fm)
  for (iz in seq_len(ps)) {
    fz <- fm
    fm <- matrix(0, nrow = nc, ncol = nc)
    for (iw in 1:ns) {
      row_idx <- ((iw - 1L) * nb + 1L):(iw * nb)
      col_idx <- ((iw - 1L) * na + 1L):(iw * na)
      fm[row_idx, col_idx] <- fz[1:nb, col_idx, drop = FALSE]
    }
    if (ns >= 2L) {
      for (ib in 2:ns) {
        col_idx <- ((ib - 1L) * nb + 1L):(ib * nb)
        fm[, col_idx] <- fm[, 1:nb, drop = FALSE]
      }
    }
    na <- na * ns; nb <- nb * ns; nc <- nc * ns
  }
  fm
}

# ---- ofn: verosimilitud, filtro (y suavizador si ks==2) ----
# Devuelve list con nll, pm, fm, skif (y skis si ks==2)
ofn <- function(th, y, ns, ps, pphi, isig, ipm,
                a = 0, b = 0, c = 0, kc = 1, ks = 1,
                izz = 1, hp) {
  # Derivados
  n <- ns^(ps + 1L)
  nk <- pphi + 1L
  capt <- length(y)
  captst <- capt - nk + 1L
  
  # --- desempaque de parámetros ---
  th <- as.numeric(th); k <- 0L
  mu  <- th[seq_len(ns)]; k <- k + ns
  phi <- if (pphi == 0L) 1 else c(1, -th[(k + 1L):(k + pphi)]); k <- k + pphi
  
  if (isig == 1L) {
    sig <- rep(th[k + 1L], n); k <- k + 1L
  } else if (isig == ns) {
    sig_prim <- th[(k + 1L):(k + ns)]; k <- k + ns
    sig <- as.numeric(t(hp[1:ns, , drop = FALSE]) %*% sig_prim)
  } else stop("isig debe ser 1 o ns.")
  if (izz == 1L) sig <- sig^2
  
  pm <- matpm(th[(k + 1L):length(th)], ns, ipm, izz)
  fm <- matf(pm, ns, ps)
  
  # --- constante dependiente del régimen ----
  stopifnot(length(phi) == ps + 1L)
  const <- as.numeric(t(kronecker(phi, mu)) %*% hp)  # 1 x n
  
  # --- construir residuales AR y densidades ---
  if (pphi == 0L) {
    Eta <- matrix(y[nk:capt], ncol = 1L)
  } else {
    Eta <- embed(y, pphi + 1L)
  }
  res   <- Eta %*% matrix(phi, ncol = 1L)              # T_eff x 1
  resM  <- matrix(res, nrow = captst, ncol = n)
  constM <- matrix(const, nrow = captst, ncol = n, byrow = TRUE)
  sigM  <- matrix(sig,   nrow = captst, ncol = n, byrow = TRUE)
  eta   <- pmin((resM - constM)^2 / sigM, 800)         # clamp a 800
  eta   <- sweep(exp(-eta / 2), 2, 1 / sqrt(sig), "*") # T_eff x n
  
  # --- ergódica inicial (LS) ---
  ap <- rbind(diag(n) - fm, matrix(1, 1, n))
  inv_tmp <- chol2inv(chol(t(ap) %*% ap))
  chsi <- pmax(rowSums(inv_tmp), 0)
  
  # --- filtro ---
  f <- 0; skif <- matrix(NA_real_, captst, n)
  for (it in 1:captst) {
    fx  <- chsi * eta[it, ]
    fit <- sum(fx); if (!is.finite(fit) || fit <= 0) fit <- 1e-300
    skif[it, ] <- fx / fit
    f <- f + log(fit)
    chsi <- as.numeric(fm %*% (fx / fit))
  }
  
  # --- ajuste bayesiano opcional ---
  if (a > 0) {
    sig_for_bayes <- sig[seq_len(ns)]
    fj <- 0
    for (ij in 1:ns)
      fj <- fj + a * log(sig_for_bayes[ij]) + b / sig_for_bayes[ij] +
      c * (mu[ij]^2 / sig_for_bayes[ij])
    f <- f - fj / 2
  }
  
  out <- list(nll = -f, pm = pm, fm = fm, skif = skif)
  
  # --- suavizador de Kim ---
  if (ks == 2L) {
    skis <- skif
    for (it in 1:(captst - 1L)) {
      t <- captst - it
      if (min(skif[t, ]) > 1e-150) {
        denom <- drop(skif[t, ] %*% t(fm)); denom[denom < 1e-150] <- 1e-150
        hk <- drop((skis[t + 1L, ] / denom) %*% fm)
        skis[t, ] <- skif[t, ] * hk
      } else {
        hk <- drop(skif[t, ] %*% t(fm))
        beta <- ifelse(hk > 1e-150, skis[t + 1L, ] / hk, 0)
        skis[t, ] <- skif[t, ] * drop(beta %*% fm)
      }
    }
    out$skis <- skis
  }
  
  if (kc == 2L) cat("\nLog likelihood:", f, "\n")
  out
}

# ---- pmth: convierte bloque de transición de x a probabilidades p(i,j) ----
pmth <- function(x, ns, ipm, ncount, nth = length(x)) {
  stopifnot(length(x) >= nth)
  if (ipm == 1L) {
    stopifnot(ns == 2L)
    x[ncount + 1L] <- x[ncount + 1L]^2 / (1 + x[ncount + 1L]^2)
    x[ncount + 2L] <- x[ncount + 2L]^2 / (1 + x[ncount + 2L]^2)
  } else if (ipm == 2L) {
    x_block <- x[(ncount + 1L):nth]
    stopifnot(length(x_block) == (ns - 1L) * ns)
    pm <- matrix(0, nrow = ns, ncol = ns)
    pm[1:(ns - 1L), ] <- matrix(x_block, nrow = ns - 1L, ncol = ns, byrow = FALSE)
    pm[ns, ] <- 1
    pm <- pm^2
    pm <- sweep(pm, 2, colSums(pm), "/")
    x[(ncount + 1L):nth] <- as.vector(pm[1:(ns - 1L), ])
  } else if (ipm == 3L) {
    idx <- (ncount + 1L):nth
    x[idx] <- x[idx]^2 / (1 + x[idx]^2)
  } else stop("ipm no soportado en pmth().")
  x
}

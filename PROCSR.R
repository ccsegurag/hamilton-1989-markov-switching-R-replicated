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
      pm <- sweep(pm, 2, colSums(pm), "/")  # normaliza **columnas** (como GAUSS)
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
    pm[3, ] <- c(1 - p11, 0,    1 - p23)  # estructura del ejemplo 5
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
  
  # --- desempaque de parámetros (orden GAUSS) ---
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
  
  # --- construir residuales AR y densidades (sin 2π) ---
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
  
  # --- suavizador de Kim (si se pide) ---
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


# ==================== FIN PROCS.R (no código suelto) ====================

# ========= Helpers TVTP =========
invlogit <- function(x) 1/(1+exp(-x))

# Construye P_t(Z) para ns=2 con logits sobre probabilidades de permanencia
# par_tvtp = c(a11, b11, a22, b22) si hay 1 exógena Z (vector)
Pt_2state <- function(z, par_tvtp) {
  a11 <- par_tvtp[1]; b11 <- par_tvtp[2]
  a22 <- par_tvtp[3]; b22 <- par_tvtp[4]
  p11 <- invlogit(a11 + b11 * z)
  p22 <- invlogit(a22 + b22 * z)
  # columnas: destino i, filas: origen j  (convención Hamilton/RATS)
  #           [ p11    1-p22
  #             1-p11   p22 ]
  matrix(c(p11, 1-p11, 1-p22, p22), nrow = 2, byrow = FALSE)
}

# --------- ofn_tvtp: NLL + filtro (+suavizador) con TVTP (ns=2) ----------
# th = (mu1, mu2, phi_1..phi_p,  sig, a11, b11, a22, b22)  [isig=1]
ofn_tvtp <- function(th, y, Z, pphi,
                     a=0, b=0, c=0,     # prior (opcional; igual que antes)
                     kc=1, ks=1) {
  
  stopifnot(pphi >= 0L, length(y) == length(Z))
  capt <- length(y)
  nk   <- pphi + 1L
  Teff <- capt - nk + 1L
  
  # ---- desempaque parámetros ----
  th <- as.numeric(th); k <- 0L
  mu  <- th[1:2];         k <- k + 2L
  phi <- if (pphi == 0L) numeric(0) else th[(k+1L):(k+pphi)]
  k   <- k + pphi
  sig <- th[k+1L]^2;      k <- k + 1L     # varianza común, izz=1→cuadrado
  par_tvtp <- th[(k+1L):(k+4L)];          # a11,b11,a22,b22
  
  # ---- construir regresores AR y residuo común (sin media) ----
  # y_t - sum phi y_{t-i}
  if (pphi > 0L) {
    Ymat <- embed(y, pphi+1L)   # última columna = y_{t-pphi}
    yt   <- Ymat[,1]
    Xphi <- Ymat[,-1,drop=FALSE]
    arpart <- drop(Xphi %*% phi)
  } else {
    yt <- y[nk:capt]; arpart <- 0
  }
  base_res <- yt - arpart                # común a todos los regímenes
  
  # ---- alinear Z_{t-1} con y_t efectivo ----
  Zlag <- Z[(nk-1L):(capt-1L)]           # misma longitud que Teff
  
  # ---- filtro con P_t(Z_{t-1}) ----
  # inicialización: ergódica de P_1 (o 0.5/0.5; aquí ergódica)
  P1   <- Pt_2state(Zlag[1], par_tvtp)
  # distribución estacionaria de P1: solve((I-P1)'; ones) en 2x2 → solución cerrada
  stat1 <- c(P1[2,2], P1[1,1]) / (P1[1,1] + P1[2,2])
  chsi  <- pmax(stat1, 1e-12); chsi <- chsi/sum(chsi)
  
  nll  <- 0
  skif <- matrix(NA_real_, Teff, 2L)
  
  for (t in 1:Teff) {
    # densidades por estado (media depende del estado)
    dens <- c( dnorm(base_res[t] - mu[1], sd = sqrt(sig)),
               dnorm(base_res[t] - mu[2], sd = sqrt(sig)) )
    # prior con TVTP
    Pt <- Pt_2state(Zlag[t], par_tvtp)
    prior <- as.numeric(Pt %*% chsi)
    # actualización
    num <- prior * dens
    fit <- sum(num); if (!is.finite(fit) || fit <= 0) fit <- 1e-300
    post <- num / fit
    skif[t,] <- post
    nll <- nll - log(fit)
    chsi <- post
  }
  
  out <- list(nll = nll, skif = skif)
  
  # ---- (opcional) suavizador con {P_{t+1}} ----
  if (ks == 2L) {
    skis <- skif
    # necesitamos también los priors del paso t+1: ξ_{t+1}^- = P_{t+1} * skif_t
    priors_next <- matrix(NA_real_, Teff-1L, 2L)
    for (t in 1:(Teff-1L)) {
      Pt1 <- Pt_2state(Zlag[t+1L], par_tvtp)
      priors_next[t,] <- as.numeric(Pt1 %*% skif[t,])
    }
    for (it in 1:(Teff-1L)) {
      t <- Teff - it
      Pt1 <- Pt_2state(Zlag[t+1L], par_tvtp)
      denom <- pmax(priors_next[t,], 1e-150)
      # factor de rescalado por componente j: skis_{t+1}(j)/denom(j)
      factor <- skis[t+1L,] / denom
      hk <- as.numeric(factor %*% Pt1)
      skis[t,] <- skif[t,] * hk
    }
    out$skis <- skis
  }
  
  if (kc == 2L) cat("\nLog-likelihood (TVTP):", -out$nll, "\n")
  out
}

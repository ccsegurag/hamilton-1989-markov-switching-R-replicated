# =================== MAXSEEK.R (versión limpia) ===================

# (opcional) fijar directorio de trabajo
setwd("C:/Users/Afsan/Desktop/Markov Switching Models - Hamilton")

# Registrar salida como en GAUSS
sink("junk")  # (acordate de sink(NULL) al final)

## ========== EXAMPLE 1 - Paper seminal de Hamilton 1989 ==========
capt <- 135L
bbbb <- scan("GNP82.DAT", quiet = TRUE)
stopifnot(length(bbbb) >= capt + 1L)
bbbb <- cbind(bbbb[1:capt], bbbb[2:(capt + 1L)])
y <- 100 * log(bbbb[, 2] / bbbb[, 1])

ns   <- 2L
ps   <- 4L
pphi <- 4L
isig <- 1L
ipm  <- 1L
a <- 0; b <- 0; c <- 0

nth <- 9L
mu1 <- 1;   mu2 <- 0
phi1 <- 0.1; phi2 <- 0; phi3 <- 0; phi4 <- 0
sig0 <- 1
p11  <- 1.5; p22 <- 1.5
th <- c(mu1, mu2, phi1, phi2, phi3, phi4, sig0, p11, p22)

# ========= pattern1 (como en GAUSS) =========
pattern1 <- function(ns, ps) {
  n  <- ns^(ps + 1L)
  na <- n / ns
  ix <- kronecker(diag(ns), matrix(1, 1, na))
  if (ps >= 1L) {
    for (i1 in 1:ps) {
      na <- na / ns
      iq <- kronecker(matrix(1, 1, ns^i1),
                      kronecker(diag(ns), matrix(1, 1, na)))
      ix <- rbind(iq, ix)
    }
  }
  ix
}
hp <- pattern1(ns, ps)

# ========= incluye PROCS (SOLO funciones) =========
source("PROCSR.R")

# ----------- bloque “startval … kc, ks …” (como GAUSS) -----------
capt <- length(y)
startval <- function() th
nk  <- pphi + 1L
izz <- 1L
n   <- ns^(ps + 1L)
kc  <- 1L
ks  <- 1L
captst <- capt - nk + 1L
id <- diag(ns)

cat("Bayesian prior used\n")
cat("a=", a, " b=", b, " c=", c, "\n", sep = "")

# ----------- objetivo (NLL) en escala de optimización -----------
obj_fn <- function(par) {
  ofn(par, y, ns, ps, pphi, isig, ipm,
      a, b, c, kc = 1, ks = 1, izz = 1, hp)$nll
}

# ----------- optimización (BFGS) -----------
x0 <- startval()
btol <- 1e-6; miter <- 150
res <- optim(par = x0, fn = obj_fn, method = "BFGS",
             control = list(maxit = miter, reltol = btol))

x_opt <- res$par          # parámetros en escala de optimización (izz=1)
f_nll <- res$value        # NLL en el óptimo

# Gradiente numérico (en escala de optimización) para imprimir como GAUSS
if (!requireNamespace("numDeriv", quietly = TRUE)) {
  install.packages("numDeriv")
}
g_opt <- numDeriv::grad(obj_fn, x_opt)

cat("\n\nMLE as parameterized for numerical optimization \n")
cat("Coefficients:\n"); print(t(x_opt))
cat("\nValue of log likelihood: ", -f_nll, "\n", sep = "")
cat("\nGradient vector:\n"); print(t(g_opt))

# ----------- Reparametrización para reporte (izz = 2) -----------
izz <- 2L
x_rep <- x_opt
ncount <- ns + pphi
x_rep[(ncount + 1L):(ncount + isig)] <- x_rep[(ncount + 1L):(ncount + isig)]^2
ncount <- ncount + isig
x_rep <- pmth(x_rep, ns = ns, ipm = ipm, ncount = ncount)  # pone p_ij en escala reporte

cat("\nVector is reparameterized to report final results as follows\n")
cat("Means for each state:\n"); print(t(x_rep[1:ns]))
if (pphi > 0L) {
  cat("Autoregressive coefficients:\n")
  print(t(x_rep[(ns + 1L):(ns + pphi)]))
}
cat("Variances:\n")
print(x_rep[(ns + pphi + 1L):(ns + pphi + isig)])

# ----------- Hessiano (escala de reporte) y SE -----------
obj_fn_report <- function(par) {
  ofn(par, y, ns, ps, pphi, isig, ipm,
      a, b, c, kc = 0, ks = 1, izz = 2, hp)$nll
}
h <- numDeriv::hessian(obj_fn_report, x_rep)
eigvals <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
if (min(eigvals) <= 0) {
  cat("Negative of Hessian is not positive definite\n")
  cat("Either not a local max, or estimates are at a boundary.\n")
} else {
  cov <- chol2inv(chol(h))
  se  <- sqrt(diag(cov))
  cat("For vector of coefficients parameterized as follows,\n"); print(t(x_rep))
  cat("the standard errors are\n"); print(t(se))
}

cat("\n-------------------------------\n")
cat("Probabilities for primitive states\n")
cat("filtered probabilities\n")

# ----------- Filtro + suavizado finales (escala reporte) -----------
kc <- 2L; ks <- 2L
out <- ofn(x_rep, y, ns, ps, pphi, isig, ipm, a, b, c, kc, ks, izz = 2, hp)

# Encabezado filtradas (ns-1 por bloque, como GAUSS)
hdr <- c("Obs")
for (tlag in 0:ps) for (i in 1:(ns - 1)) hdr <- c(hdr, sprintf("P(st-%d=%d)", tlag, i))
cat(paste(hdr, collapse = " "), "\n")

E <- kronecker(diag(ps + 1L), diag(ns)[, 1:(ns - 1L), drop = FALSE])
skif_prim <- (out$skif %*% t(hp)) %*% E
skif_print <- cbind(seq.int(from = nk, by = 1L, length.out = captst),
                    round(skif_prim, 4))
print(skif_print, row.names = FALSE)

cat("\nsmoothed probabilities\n")
hdr2 <- c("Obs", paste0("P(st=", 1:ns, ")"))
cat(paste(hdr2, collapse = " "), "\n")

skis_prim <- out$skis %*% t(hp)              # tomar bloque actual s_t
smoothed_print <- cbind(seq.int(from = nk, by = 1L, length.out = captst),
                        round(skis_prim[, 1:ns, drop = TRUE], 4))
print(smoothed_print, row.names = FALSE)

# Cerrar archivo de salida
sink(NULL)

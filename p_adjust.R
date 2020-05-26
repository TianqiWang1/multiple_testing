p.adjust.methods <-
  c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    "ABH", "MABH", "storey", "TST", "MTST",  # <<<--- the new methods
    "none")

p.adjust <- function(p, method = p.adjust.methods, n = length(p),
                     alpha = 0.05, lambda = 0.5)
{
  ## Methods 'Hommel', 'BH', 'BY' and speed improvements
  ## contributed by Gordon Smyth
  method <- match.arg(method)
  if(method == "fdr") method <- "BH"	# back compatibility
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if(all(nna <- !is.na(p))) nna <- TRUE
  p <- p[nna]
  lp <- length(p)
  stopifnot(n >= lp)
  if (n <= 1) return(p0)
  if (n == 2 && method == "hommel") method <- "hochberg"
  
  p0[nna] <-
    switch(method,
           bonferroni = pmin(1, n * p),
           holm = {
             i <- seq_len(lp)
             o <- order(p)
             ro <- order(o)
             pmin(1, cummax( (n - i + 1L) * p[o] ))[ro]
           },
           hommel = { ## needs n-1 >= 2 in for() below
             if(n > lp) p <- c(p, rep.int(1, n-lp))
             i <- seq_len(n)
             o <- order(p)
             p <- p[o]
             ro <- order(o)
             q <- pa <- rep.int( min(n*p/i), n)
             for (j in (n-1):2) {
               ij <- seq_len(n-j+1)
               i2 <- (n-j+2):n
               q1 <- min(j*p[i2]/(2:j))
               q[ij] <- pmin(j*p[ij], q1)
               q[i2] <- q[n-j+1]
               pa <- pmax(pa,q)
             }
             pmax(pa,p)[if(lp < n) ro[1:lp] else ro]
             
           },
           hochberg = {
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             pmin(1, cummin( (n - i + 1L) * p[o] ))[ro]
             
           },
           BH = {
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             pmin(1, cummin( n / i * p[o] ))[ro]
           },
           
           BY = {
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             q <- sum(1L/(1L:n))
             pmin(1, cummin(q * n / i * p[o]))[ro]
             
           },
           
           ABH = {
             o <- order(p)
             i <- 1L:lp
             h0.m <- (n+1-i)/(1-p[o])
             idx <- which(diff(h0.m)>0)
             if (length(idx) != 0) {
               m_prev <- h0.m[min(idx) + 1]
               m <- ceiling(min(m_prev, lp))
             }
             else {
               m <- lp
             }
             
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             pmin(1, cummin( m / i * p[o] ))[ro]
           },
           
           TST = {
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             p_adj <- pmin(1, cummin( n / i * p[o] ))[ro]
             m <- sum(p_adj >= alpha/(1+alpha),na.rm=TRUE)
             if (m != 0 & m != n) {
               p_adj <- pmin(1, cummin( m / i * p[o] ))[ro]
             }
             p_adj
           },
           
           MABH = {
             m <- (n/2)/(1-median(p))
             m <- ceiling(min(m, lp))
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             pmin(1, cummin( m / i * p[o] ))[ro]
           },
           
           storey = {
             r <- sum(p <= lambda, na.rm=TRUE)
             m <- (n - r + 1)/(1 - lambda)
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             pmin(1, cummin( m / i * p[o] ))[ro]
           },
           
           MTST = {
             i <- lp:1L
             o <- order(p, decreasing = TRUE)
             ro <- order(o)
             p_adj <- pmin(1, cummin( n / i * p[o] ))[ro]
             m <- length(which(p_adj >= alpha))
             if (m != 0 & m != n) {
               p_adj <- pmin(1, cummin( m / i * p[o] * (1+alpha)))[ro]
             }
             p_adj
           },
           none = p)
  p0
}
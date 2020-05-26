library(multtest)
library(sfsmisc)

alpha <- 0.05
n <- 1000

############
# This function outputs FDR, type 1 error rate, and power,
# given the adjusted_p_value, m, n, and true_frac as input
# It's called inside the sim() function
############

compute_results<- function(adjusted_p,m, n, true_frac, ret){ 
  adjusted_p <- matrix(adjusted_p, n, m)
  # rejection index
  rejected_idx <- apply(adjusted_p, 1, function(x) which(x < alpha))
  # number of false rejection
  num_false_rejected <- lapply(rejected_idx, function(x) length(which(x <= true_frac*m)))
  # number of rejections
  num_rejected <- lapply(rejected_idx, length)
  
  res <- c()
  
  # FDR is #false rejection/#rejection
  if ("FDR" %in% ret) {
    num_rejected_adj <- lapply(num_rejected, function(x) max(x,1))
    fdr <- Map('/', num_false_rejected, num_rejected_adj)
    FDR <- round(mean(unlist(fdr)),3)
    res <- c(res, FDR)
  }
  
  # Type 1 error rate is #false rejection/#true null
  if ("Type_1" %in% ret) {
    type_1 <- lapply(num_false_rejected, function(x) x/(true_frac*m))
    type_1_error_rate <- round(mean(unlist(type_1)), 3)
    res <- c(res, type_1_error_rate)
  }
  
  # power is #rejected false null/#false null
  if ("Power" %in% ret) {
    num_true_rejected <- Map('-', num_rejected, num_false_rejected)
    pow <- lapply(num_true_rejected, function(x) x/((1-true_frac)*m))
    powers <- round(mean(unlist(pow)),3)
    res <- c(res, powers)
  }
  return(res)
}

##################
# This function takes in the simulation result
# and plots result vs. m
##################

plot_single <- function(sim_res, m_vals, cor, true_frac, procs = procs, cex=0.4, legend=TRUE, ret="fdr", ymin, ymax) {
  colors <- rainbow(dim(sim_res)[1])
  linetype <- c(1:dim(sim_res)[1])
  
  plot(1,xlim = c(1, dim(sim_res)[2]), ylim = c(ymin, ymax), type='n', xlab = 'm', ylab= ret,xaxt = "n" )
  
  for (i in 1:dim(sim_res)[1]) {
    lines(c(1:dim(sim_res)[2]), sim_res[i,], col=colors[i], lty=linetype[i])
  }
  axis(1, at=c(1:dim(sim_res)[2]), labels=as.character(m_vals)) 
  
  if (legend) {
    legend("bottomright", legend=c(procs,"oracle"),cex=cex, col=colors,lty=linetype,text.font=4, bty="n")
    #legend("topleft",legend=c(procs,"oracle"), cex=cex,col=colors, lty=linetype, text.font=4, inset=c(1,0), xpd=TRUE, bty="n")
  }
  title(main=sprintf("cor = %s, m0/m = %s", cor, true_frac), font.main=1)
}

##########
# This function outputs the result for selected procedures as a matrix,
# given m, cor, true_frac, and n as input
# We could save the matrix as .csv or .Rdata
##########

sim <- function(cor, true_frac, n, m, procs = c("bonferroni", "holm","hochberg", "BH", "BY","ABH", "MABH","storey", "TST", "MTST"), ret = c("FDR", "Type_1", "Power"), mu_config = 1, save=FALSE) {
  set.seed(123)
  # generate Z
  Z0 <- matrix(rep(rnorm(n,0,1),m),n,m)
  Z <- matrix(rnorm(n*m, 0, 1), n, m)
  true_len <- as.integer(m*true_frac)
  
  if (mu_config == 1) {
    # "all at 5" configuration
    mu <- matrix(c(rep(0, true_len*n), rep(5, (m-true_len)*n)), n, m)
  }  else {
    # "1 2 3 4" configuration
    left_mu <- matrix(c(rep(0, true_len*n)), n, true_len)
    right_mu <- t(matrix(rep(c(1,2,3,4), (m-true_len)*n/4), m-true_len,n))
    mu <- cbind(left_mu, right_mu)
  }
  # generate Y
  Y <- sqrt(cor)*Z0 + sqrt(1-cor)*Z + mu
  # generate P
  P <- 2*(1-pnorm(abs(Y)))
  # matrix for storing the result
  mat <- matrix(0,0,length(ret))
  # for each procedure, compute the result
  for(proc in procs) {
    p_vals <- p.adjust(P, method=proc)
    cur_res <- compute_results(p_vals, m, n, true_frac, ret)
    mat <- rbind(mat, cur_res)
  }
  # compute the result for oracle
  p_adj_oracle <- p.adjust(P*true_frac, method="BH")
  oracle_res <- compute_results(p_adj_oracle, m, n, true_frac, ret)
  mat <- rbind(mat, oracle_res)
  rownames(mat) <- c(procs,"oracle")
  colnames(mat) <- ret
  
  if (save) {
    save_csv = sprintf("mat_m%s_cor%s_frac%s.csv", m, cor, true_frac)
    write.csv(mat, file = save_csv)
  }
  return(mat)
}


################
# This function takes in a list of cor values, m values, and true fraction values
# the result is stored in a single matrix
# can store the result as .csv or .Rdata
# can plot result vs. m
###############

mult_sim <- function(cor_vals, true_frac_vals, m_vals,ret = c("FDR", "Type_1", "Power"), procs = c("bonferroni", "hochberg", "holm", "BH", "BY","ABH", "MABH","storey", "TST", "MTST") ,mu_config = 1, save=FALSE, plot=FALSE, ymin=0, ymax=0.06, cex=0.4) {
  
  mat <- matrix(0,length(procs)+1,0)
  
  for (cor in cor_vals) {
    
    for (true_frac in true_frac_vals) {
      
      plot_mat <- matrix(0, length(procs)+1, 0)
      
      for (m in m_vals) {
        
        res <- sim(cor,true_frac,n, m, ret =ret, procs = procs, mu_config=mu_config)
        plot_mat <- cbind(plot_mat, res)
        mat <- cbind(mat, res)
        
      }
      
      if (plot) {
        lg = true_frac == true_frac_vals[length(true_frac_vals)]
        plot_single(plot_mat, m_vals, cor, true_frac, procs, ret=ret, legend=lg, ymin=ymin, ymax = ymax, cex=cex)
      }
      
    }
  }
  if (save) {
    save_csv = sprintf("table_1.csv", cor, true_frac)
    write.csv(mat, file = save_csv)
  }
  return(mat)
} 

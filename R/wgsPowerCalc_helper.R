#' Tests the power to detect a single locus using case control data varying relative risk.
#'
#' @param f_gene Number of possible functional rare variants per gene (numeric).
#' @param s Number of possible rare variants per gene (numeric).
#' @param AF_bar Average minor allele frequency (numeric, 0-1).
#' @param N_rep Number of replicate simulations (integer).
#' @param R Maximum relative risk (numeric).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param K Prevalence of the disorder (numeric, 0-1).
#' @param N Sample size (integer).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_cc_locus_single p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the relative risk is varied.
#' @examples
#' f_gene <- 123
#' s <- 615
#' AF_bar <- 0.001
#' N_rep <- 50
#' R <- 2.5
#' f <- 0.2
#' K <- 0.01
#' N <- 20000
#' r <- 1
#' p_thres_cc_locus_single <- 0.05 / 3000000000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotCcLocusSingleByRelativeRisk()
#' @export
plotCcLocusSingleByRelativeRisk <- function(f_gene, s, AF_bar, N_rep, R, f, K, N, r, p_thres_cc_locus_single=0.05, name="Test", col="#e41a1c"){

  if(f_gene > s){
    stop("Number of functonal variant must be less than total number of variants.", call. = FALSE)
  }

  # Creates an exponential distribution of allele frequencies for each possible rare variant by each replicate
  AF_raw <- rexp(s*N_rep, 1/AF_bar)
  AF_raw <- matrix(AF_raw, s, N_rep)

  # Samples the exponential distribution to get allele frequencies for each possible rare variant by each replicate
  AF_save <- lapply(c(1:N_rep), function(i_rep, AF_raw){
    prob <- 1/AF_raw[,i_rep];
    prob <- prob/sum(prob);
    ret <- sample(AF_raw[,i_rep], prob=prob)
    ret}, AF_raw)
  AF_save <- do.call("cbind", AF_save)

  # Relative risk values to be tested
  rr <- seq(1, R, by=0.05)

  # For each replicate estimate power at each relative risk
  power_single <- c()
  for (i_rep in c(1:N_rep)){
    power_single <- cbind(power_single, unlist(lapply(rr, function(R, f, AF, K, N, r, p_thres_cc_locus_single) {getSingleVarLocusPower(f, AF, R, K, N, r, p_thres_cc_locus_single)}, f_gene, AF_save[,i_rep],  K, N, r, p_thres_cc_locus_single)))
  }

  # Print result to a PDF
  nnName <- paste(name,"_CC_LocusSingle_RelRisk_s_",s,"_f_",f,".pdf", sep = "")
  pdf(nnName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Relative risk",  ylab="Power", xlim=range(rr), ylim=c(0,1))
  lines(rr, apply(power_single, 1, mean), lwd=3, col=col)
  dev.off()
}

#' Tests the power to detect a single locus using case control data varying sample size
#'
#' @param f_gene Number of possible functional rare variants per gene (numeric).
#' @param s Number of possible rare variants per gene (numeric).
#' @param AF_bar Average minor allele frequency (numeric, 0-1).
#' @param N_rep Number of replicate simulations (integer).
#' @param R Relative risk (numeric).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param K Prevalence of the disorder (numeric, 0-1).
#' @param N Maximum sample size (integer).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_cc_locus_single p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the sample size is varied.
#' @examples
#' f_gene <- 123
#' s <- 615
#' AF_bar <- 0.001
#' N_rep <- 50
#' R <- 2.5
#' f <- 0.2
#' K <- 0.01
#' N <- 20000
#' r <- 1
#' p_thres_cc_locus_single <- 0.05 / 3000000000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotCcLocusSingleBySampleSize()
#' @export
plotCcLocusSingleBySampleSize <- function(f_gene, s, AF_bar, N_rep, R, f, K, N, r, p_thres_cc_locus_single=0.05, name="Test", col="#e41a1c"){

  if(f_gene > s){
    stop("Number of functonal variant must be less than total number of variants.", call. = FALSE)
  }

  # Creates an exponential distribution of allele frequencies for each possible rare variant by each replicate
  AF_raw <- rexp(s*N_rep, 1/AF_bar)
  AF_raw <- matrix(AF_raw, s, N_rep)

  # Samples the exponential distribution to get allele frequencies for each possible rare variant by each replicate
  AF_save <- lapply(c(1:N_rep), function(i_rep, AF_raw){
    prob <- 1/AF_raw[,i_rep];
    prob <- prob/sum(prob);
    ret <- sample(AF_raw[,i_rep], prob=prob)
    ret}, AF_raw)
  AF_save <- do.call("cbind", AF_save)

  # Sample size values to be tested
  nn <- seq(100, N, by=10000)

  # For each replicate estimate power at each sample size
  power_single <- c()
  for (i_rep in c(1:N_rep)){
    power_single <- cbind(power_single, unlist(lapply(nn, function(N, f, AF, K, R, r, p_thres_cc_locus_single) {getSingleVarLocusPower(f, AF, R, K, N, r, p_thres_cc_locus_single)}, f_gene, AF_save[,i_rep],  K, R, r, p_thres_cc_locus_single)))
  }

  # Print result to a PDF
  nnName <- paste(name,"_CC_LocusSingle_SampSize_s_",s,"_f_",f,".pdf", sep = "")
  pdf(nnName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Sample Size",  ylab="Power", xlim=range(nn), ylim=c(0,1))
  lines(nn, apply(power_single, 1, mean), lwd=3, col=col)
  dev.off()
}

#' Tests the power to detect a gene locus (i.e. locus defined by a gene rather than variant) using case control data varying relative risk.
#'
#' @param f_gene Number of possible functional rare variants per gene (numeric).
#' @param s Number of possible rare variants per gene (numeric).
#' @param AF_bar Average minor allele frequency (numeric, 0-1).
#' @param N_rep Number of replicate simulations (integer).
#' @param R Maximum relative risk (numeric).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param K Prevalence of the disorder (numeric, 0-1).
#' @param N Sample size (integer).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_cc_locus_multi p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the relative risk is varied.
#' @examples
#' f_gene <- 123
#' s <- 615
#' AF_bar <- 0.001
#' N_rep <- 50
#' R <- 2.5
#' f <- 0.2
#' K <- 0.01
#' N <- 20000
#' r <- 1
#' p_thres_cc_locus_multi <- 0.05 / 20000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotCcLocusMultiByRelativeRisk()
#' @export
plotCcLocusMultiByRelativeRisk <- function(f_gene, s, AF_bar, N_rep, R, f, K, N, r, p_thres_cc_locus_multi=0.05, name="Test", col="#e41a1c"){

  if(f_gene > s){
    stop("Number of functonal variant must be less than total number of variants.", call. = FALSE)
  }

  # Creates an exponential distribution of allele frequencies for each possible rare variant by each replicate
  AF_raw <- rexp(s*N_rep, 1/AF_bar)
  AF_raw <- matrix(AF_raw, s, N_rep)

  # Samples the exponential distribution to get allele frequencies for each possible rare variant by each replicate
  AF_save <- lapply(c(1:N_rep), function(i_rep, AF_raw){
    prob <- 1/AF_raw[,i_rep];
    prob <- prob/sum(prob);
    ret <- sample(AF_raw[,i_rep], prob=prob)
    ret}, AF_raw)
  AF_save <- do.call("cbind", AF_save)

  # Calculates the sum of all heterozygous allele freuqencies (p * 1-p) for each replicate
  sum_var <- apply(AF_save, 2, function(x){ sum(x*(1-x) ) })

  # Relative risk values to be tested
  rr <- seq(1, R, by=0.05)

  # For each replicate estimate power at each relative risk
  power_multi <- c()
  for (i_rep in c(1:N_rep)){
    power_multi <- cbind(power_multi, unlist(lapply(rr, function(R, f, AF, sum_var, K, N, r, p_thres_cc_locus_multi){getMultiVarLocusPower(f, AF, sum_var, R, K, N, r, p_thres_cc_locus_multi)}, f_gene, AF_save[,i_rep], sum_var[i_rep], K, N, r, p_thres_cc_locus_multi )))
  }

  # Print result to a PDF
  nnName <- paste(name,"_CC_LocusMulti_RelRisk_s_",s,"_f_",f,".pdf", sep = "")
  pdf(nnName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Relative risk",  ylab="Power", xlim=range(rr), ylim=c(0,1))
  lines(rr, apply(power_multi, 1, mean), lwd=3, col=col)
  dev.off()
}

#' Tests the power to detect a gene locus (i.e. locus defined by a gene rather than variant) using case control data varying sample size.
#'
#' @param f_gene Number of possible functional rare variants per gene (numeric).
#' @param s Number of possible rare variants per gene (numeric).
#' @param AF_bar Average minor allele frequency (numeric, 0-1).
#' @param N_rep Number of replicate simulations (integer).
#' @param R Relative risk (numeric).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param K Prevalence of the disorder (numeric, 0-1).
#' @param N Maximum sample size (integer).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_cc_locus_multi p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the sample size is varied.
#' @examples
#' f_gene <- 123
#' s <- 615
#' AF_bar <- 0.001
#' N_rep <- 50
#' R <- 2.5
#' f <- 0.2
#' K <- 0.01
#' N <- 20000
#' r <- 1
#' p_thres_cc_locus_multi <- 0.05 / 20000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotCcLocusMultiBySampleSize()
#' @export
plotCcLocusMultiBySampleSize <- function(f_gene, s, AF_bar, N_rep, R, f, K, N, r, p_thres_cc_locus_multi=0.05, name="Test", col="#e41a1c"){

  if(f_gene > s){
    stop("Number of functonal variant must be less than total number of variants.", call. = FALSE)
  }

  # Creates an exponential distribution of allele frequencies for each possible rare variant by each replicate
  AF_raw <- rexp(s*N_rep, 1/AF_bar)
  AF_raw <- matrix(AF_raw, s, N_rep)

  # Samples the exponential distribution to get allele frequencies for each possible rare variant by each replicate
  AF_save <- lapply(c(1:N_rep), function(i_rep, AF_raw){
    prob <- 1/AF_raw[,i_rep];
    prob <- prob/sum(prob);
    ret <- sample(AF_raw[,i_rep], prob=prob)
    ret}, AF_raw)
  AF_save <- do.call("cbind", AF_save)

  # Calculates the sum of all heterozygous allele freuqencies (p * 1-p) for each replicate
  sum_var <- apply(AF_save, 2, function(x){ sum(x*(1-x) ) })

  # Sample size values to be tested
  nn <- seq(100, N, by=5000)

  # For each replicate estimate power at each sample size
  power_multi <- c()
  for (i_rep in c(1:N_rep)){
    power_multi <- cbind(power_multi, unlist(lapply(nn, function(N, f, AF, sum_var, K, R, r, p_thres_cc_locus_multi){getMultiVarLocusPower(f, AF, sum_var, R, K, N, r, p_thres_cc_locus_multi)}, f_gene, AF_save[,i_rep], sum_var[i_rep], K, R, r, p_thres_cc_locus_multi )))
  }

  # Print result to a PDF
  nnName <- paste(name,"_CC_LocusMulti_SampSize_s_",s,"_f_",f,".pdf", sep = "")
  pdf(nnName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Sample size",  ylab="Power", xlim=range(nn), ylim=c(0,1))
  lines(nn, apply(power_multi, 1, mean), lwd=3, col=col)
  dev.off()
}


#' Tests the power to detect a burden difference using case control data varying relative risk.
#'
#' @param R Maximum relative risk (numeric).
#' @param N Sample size (integer).
#' @param q Number of functional rare variants per person in seleceted regions (integer).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param p_thres_cc_burden p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the relative risk is varied.
#' @examples
#' R <- 2.5
#' N <- 20000
#' q <- 126
#' f <- 0.2
#' p_thres_cc_burden <- 0.05 / 1000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotCcBurdenByRelativeRisk()
#' @export
plotCcBurdenByRelativeRisk <- function (R, N, q, f, p_thres_cc_burden=0.05, name="Test", col="#e41a1c"){
  rr <- seq(1, R, by=.025)

  power_burden <- c()
  power_burden <- cbind(power_burden, unlist(lapply(rr, function(R, N, q, f){getCaseConBurdenPower(R, N, q, f)}, N, q, f)))

  # Print result to a PDF
  rrName <- paste(name,"_CC_Burden_RelRisk_q_",q,"_f_",f,".pdf", sep = "")
  pdf(rrName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Relative risk", ylab="Power", xlim=range(rr), ylim=c(0,1))
  lines(rr, power_burden, lwd=3, col=col)
  dev.off()
}

#' Tests the power to detect a burden difference using case control data varying sample size.
#'
#' @param R Relative risk (numeric).
#' @param N Maximum sample size (integer).
#' @param q Number of functional rare variants per person in seleceted regions (integer).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param p_thres_cc_burden p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the sample size is varied.
#' @examples
#' R <- 2.5
#' N <- 20000
#' q <- 126
#' f <- 0.2
#' p_thres_cc_burden <- 0.05 / 1000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotCcBurdenBySampleSize()
#' @export
plotCcBurdenBySampleSize <- function (R, N, q, f, p_thres_cc_burden=0.05, name="Test", col="#e41a1c"){
  nn <- seq(1, N, by=1000)

  power_burden <- c()
  power_burden <- cbind(power_burden, unlist(lapply(nn, function(N, R, q, f){getCaseConBurdenPower(R, N, q, f)}, R, q, f)))

  # Print result to a PDF
  nnName <- paste(name,"_CC_Burden_SampSize_q_",q,"_f_",f,".pdf", sep = "")
  pdf(nnName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Sample size", ylab="Power", xlim=range(nn), ylim=c(0,1))
  lines(nn, power_burden, lwd=3, col=col)
  dev.off()
}

#' Tests the power to detect a burden difference using de novo data varying relative risk.
#'
#' @param R Maximum relative risk (numeric).
#' @param N Sample size (integer).
#' @param q Number of functional rare variants per person in seleceted regions (integer).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param p_thres_denovo_burden p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the relative risk is varied.
#' @examples
#' R <- 25
#' N <- 5000
#' q <- 0.0249
#' f <- 0.2
#' p_thres_denovo_burden <- 0.05 / 1000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotDnBurdenByRelativeRisk()
#' @export
plotDnBurdenByRelativeRisk <- function (R, N, q, f, p_thres_denovo_burden=0.05, name="Test", col="#e41a1c"){
  rr <- seq(1, R, by=.1)

  power_burden <- c()
  power_burden <- cbind(power_burden, unlist(lapply(rr, function(R, N, q, f){getDeNovoBurdenPower(R, N, q, f)}, N, q, f)))

  # Print result to a PDF
  rrName <- paste(name,"_DN_Burden_RelRisk_q_",q,"_f_",f,".pdf", sep = "")
  pdf(rrName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Relative risk", ylab="Power", xlim=range(rr), ylim=c(0,1))
  lines(rr, power_burden, lwd=3, col=col)
  dev.off()
}

#' Tests the power to detect a burden difference using de novo data varying sample size.
#'
#' @param R Relative risk (numeric).
#' @param N Maximum Sample size (integer).
#' @param q Number of functional rare variants per person in seleceted regions (integer).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param p_thres_denovo_burden p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the sample size is varied.
#' @examples
#' R <- 25
#' N <- 5000
#' q <- 0.0249
#' f <- 0.2
#' p_thres_denovo_burden <- 0.05 / 1000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotDnBurdenBySampleSize()
#' @export
plotDnBurdenBySampleSize <- function (R, N, q, f, p_thres_denovo_burden=0.05, name="Test", col="#e41a1c"){
  nn <- seq(1, N, by=100)

  power_burden <- c()
  power_burden <- cbind(power_burden, unlist(lapply(nn, function(N, R, q, f){getDeNovoBurdenPower(R, N, q, f)}, R, q, f)))

  # Print result to a PDF
  nnName <- paste(name,"_DN_Burden_SampSize_q_",q,"_f_",f,".pdf", sep = "")
  pdf(nnName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Sample size", ylab="Power", xlim=range(nn), ylim=c(0,1))
  lines(nn, power_burden, lwd=3, col=col)
  dev.off()
}


#' Tests the power to detect a gene locus using de novo data varying relative risk.
#'
#' @param R Maximum relative risk (numeric).
#' @param N Sample size (integer).
#' @param q Number of functional rare variants per person in seleceted regions (integer).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_denovo_locus p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the relative risk is varied.
#' @examples
#' R <- 25
#' N <- 5000
#' q <- 0.0249
#' f <- 0.2
#' r <- 1
#' p_thres_denovo_locus <- 0.05 / 20000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotDnLocusByRelativeRisk()
#' @export
plotDnLocusByRelativeRisk <- function(R, N, q, f, r, p_thres_denovo_locus=0.05, name="Test", col="#e41a1c"){
  rr <- seq(1, R, by=.2)

  power_para_selected <- getDeNovoLocusPowerParametric(rr, q, f, N, r, p_thres_denovo_locus)

  # Print result to a PDF
  nnName <- paste(name,"_DN_Locus_RelRisk_q_",q,"_f_",f,".pdf", sep = "")
  pdf(nnName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Relative risk", ylab="Power", xlim=range(rr), ylim=c(0,1))
  lines(rr, power_para_selected, lwd=3, col=col)
  dev.off()
}

#' Tests the power to detect a gene locus using de novo data varying sample size.
#'
#' @param R Relative risk (numeric).
#' @param N Maximum sample size (integer).
#' @param q Number of functional rare variants per person in seleceted regions (integer).
#' @param f Proportion of functional rare variants that mediate risk for the disorder (numeric, 0-1).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_denovo_locus p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @param name Prefix for plot filenames (text).
#' @param col Color of the line (hex).
#' @return A PDF showing the estimated statistical power as the sample size is varied.
#' @examples
#' R <- 25
#' N <- 5000
#' q <- 0.0249
#' f <- 0.2
#' r <- 1
#' p_thres_denovo_locus <- 0.05 / 20000
#' name <- "Test1"
#' col <- "#e41a1c"
#' plotDnLocusBySampleSize()
#' @export
plotDnLocusBySampleSize <- function(R, N, q, f, r, p_thres_denovo_locus=0.05, name="Test", col="#e41a1c"){
  nn <- seq(100, N, by=100)

  power_para_selected <- getDeNovoLocusPowerParametric(R, q , f, nn, r, p_thres_denovo_locus)

  # Print result to a PDF
  nnName <- paste(name,"_DN_Locus_SampSize_q_",q,"_f_",f,".pdf", sep = "")
  pdf(nnName, onefile = FALSE, paper = "special", height = 6, width = 6)
  plot(0, 0, type="n", xlab="Sample size", ylab="Power", xlim=range(nn), ylim=c(0,1))
  lines(nn, power_para_selected, lwd=3, col=col)
  dev.off()
}

#' Calculate statistical power in a case control burden analysis.
#'
#' @param R Relative Risk (numeric).
#' @param N Sample Size (integer).
#' @param q Total number of variants within selected regions per control individual (numeric).
#' @param f Percent of variants within selected regions that mediate disease risk (numeric).
#' @return The statistical power in a case control cohort.
#' @examples
#' getCaseConBurdenPower(2.5, 20000, 1.26, 0.2)
#' @export
getCaseConBurdenPower <- function (R, N, q, f){
  # Number of variants in controls
  conNum <- q * N
  # Number of variants in cases
  caseNum <- ( q * N * f * R ) + ( q * N * (1-f) )
  sd <- 0.91 * conNum/0.67 # Estimate based on SSC exome data
  burden <- caseNum/conNum
  r <- 1 #Case:Control ratio (numeric).
  
  # Estimate power
  z <- (caseNum-conNum)/(sd*sqrt((1+1/r)/N))
  power <- pnorm(z-qnorm(1-p_thres_cc_burden/2))+pnorm(-z-qnorm(1-p_thres_cc_burden/2))
  power
}

#' Calculate statistical power in a de novo burden analysis.
#'
#' @param R Relative Risk (numeric).
#' @param N Sample Size (integer).
#' @param q Total number of variants within selected regions per control individual (numeric).
#' @param f Percent of variants within selected regions that mediate disease risk (numeric).
#' @return The statistical power in a de novo cohort.
#' @examples
#' getDeNovoBurdenPower(25, 5000, 0.0249, 0.2)
#' @export
getDeNovoBurdenPower <- function (R, N, q, f){
  # Number of variants in controls
  conNum <- q * N
  # Number of variants in cases
  caseNum <- ( q * N * f * R ) + ( q * N * (1-f) )
  sd <- 0.31 * conNum/0.0966 # Estimate based on SSC exome data
  burden <- caseNum/conNum
  r <- 1 #Case:Control ratio (numeric).

  # Estimate power
  z <- (caseNum-conNum)/(sd*sqrt((1+1/r)/N))
  power <- pnorm(z-qnorm(1-p_thres_denovo_burden/2))+pnorm(-z-qnorm(1-p_thres_denovo_burden/2))
  power
}

#' Calculates the chi-square distribution with the non-centrality parameter χ2 in cases and controls.
#'
#' @param q Allele frequency (numeric, 0-1).
#' @param R Relative Risk (numeric).
#' @param K Prevalence of the disorder (numeric, 0-1).
#' @param CaseN Number of cases (integer).
#' @param ConN Number of controls (integer).
#' @return The chi-squared value (numeric).
#' @examples
#' NCPgen(1e-05, 2.5, 0.01, 10000, 10000)
#' @export
NCPgen <- function(q, R, K, CaseN, ConN) {
  p <- 1-q;

  p_d_aa <- K / (p*p + 2*p*q*R + q*q*R*R);
  p_d_ab <- p_d_aa * R;
  p_d_bb <- p_d_aa * R * R;

  p_aa_d <- p_d_aa * p * p/ K ;
  p_ab_d <- p_d_ab * p * q * 2 / K ;
  p_bb_d <- p_d_bb * q * q/ K ;

  p_aa_u <- (p*p - K * p_aa_d)  / (1-K) ;
  p_ab_u <- (2*p*q - K * p_ab_d) / (1-K) ;
  p_bb_u <- (q*q - K * p_bb_d)  / (1-K) ;

  p_case <- p_aa_d + 0.5*p_ab_d;
  p_con <- p_aa_u + 0.5*p_ab_u;
  q_case <- 1 - p_case ;
  q_con <- 1 - p_con ;

  chi <- summary(as.table(matrix(c(2*p_case*CaseN,2*p_con*ConN,2*q_case*CaseN,2*q_con*ConN),2,2)))$statistic
  chi
}

#' Calculates the power of rare variant locus discovery in cases and controls.
#'
#' @param f Proportion of variants that mediate risk (numeric, 0-1).
#' @param AF Allele frequency (numeric, 0-1).
#' @param R Relative Risk (numeric).
#' @param K Prevalence of the disorder (numeric, 0-1).
#' @param N Sample Size (integer).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_cc_locus_single p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @return The statistical power (numeric, 0-1).
#' @examples
#' getSingleVarLocusPower(0.2, 1e-4, 2.5, 0.01, 20000, 1, 1.666667e-11)
#' @export
getSingleVarLocusPower <- function(f, AF, R, K, N, r, p_thres_cc_locus_single){
  chi <- unlist(lapply(AF[1:f], function(i){ NCPgen(i, R, K, N*r/(1+r), N/(1+r)) }))
  power_snp <- pchisq(qchisq(p_thres_cc_locus_single, 1, low=F), 1, ncp=chi, low=F)
  ret <- 1-prod(1-power_snp)
  ret
}

#' Calculates the power of rare variant locus discovery across a gene in cases and controls.
#'
#' @param f Proportion of variants that mediate risk (numeric, 0-1).
#' @param AF Allele frequency (numeric, 0-1).
#' @param sum_var The sum of all heterozygous allele freuqencies (p * 1-p) for each replicates (numeric)
#' @param R Relative Risk (numeric).
#' @param K Prevalence of the disorder (numeric, 0-1).
#' @param N Sample Size (integer).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_cc_locus_multi p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @return The statistical power (numeric, 0-1).
#' @examples
#' getMultiVarLocusPower(0.2, 1e-4, 0.6, 2.5, 0.01, 20000, 1, 2.5e-06)
#' @export
getMultiVarLocusPower <- function(f, AF, sum_var, R, K, N, r, p_thres_cc_locus_multi){
  ratio <- sum(AF[1:f] * (1-AF[1:f]) ) / sum_var
  chi <- NCPgen(sum(AF), exp(log(R)*ratio), K, N*r/(1+r), N/(1+r))
  power <- pchisq(qchisq(p_thres_cc_locus_multi, 1, low=F), 1, ncp=chi , low=F)
  power
}

#' Calculates the power of de novo mutation burden in de novo cohorts.
#'
#' @param R Relative Risk (numeric).
#' @param q Total number of de novo mutations within selected regions per control individual (numeric).
#' @param f Proportion of de novo mutations that mediate risk (numeric, 0-1).
#' @param N Sample Size (integer).
#' @param r Case:Control ratio (numeric).
#' @param p_thres_denovo_locus p-value threshold after correction for multiple comparisons (numeric, 0-1).
#' @return The statistical power (numeric, 0-1).
#' @examples
#' getDeNovoLocusPowerParametric(25, 0.0249, 0.2, 5000, 1, 5e-05)
#' @export
getDeNovoLocusPowerParametric <- function(R, q, f, N, r, p_thres_denovo_burden){
  MutCases <- N*r/(1+r)
  MutCon <- N/(1+r)
  x1 <- q * 1* (1-f) +  q * R * f
  x0 <- q * 1
  chi <- (x1-x0)^2/(x1/MutCases+x0/MutCon)
  power <- pchisq(qchisq(p_thres_denovo_burden, 1, low=F), 1, ncp=chi, low=F)
  power
}


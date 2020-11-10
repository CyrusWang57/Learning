library(corrplot) 
library(effects) 
library(nlme)
setwd("E:\wangc\Documents\R")
source('toolbox_coexistence.R')
source('toolbox_figure.R')
library(codyn) 
library(MASS)
##ρ=-0.1====
c<-2
phi_e<-0.5
#var<-0.01^2
KK<-20000
LV_mono <- function(r, a, N0, steps,var_e,var_d) {
  N <- matrix(c(N0, rep(numeric(steps - 1), length(r))), steps, length(r), byrow=T)
  COV <- matrix(rep(-0.1*var_e,sr_i^2),nrow=sr_i,ncol=sr_i)
  diag(COV) <- var_e #construct a diagonal matrix.对角矩阵：一个主对角线之外的元素皆为0的矩阵
  COV2<- matrix(rep(-0.1*var_d,sr_i^2),nrow=sr_i,ncol=sr_i)
  diag(COV2) <- var_d
  for(t in 1: (steps - 1)) { 
    e<-mvrnorm(1,rep(0,length(r)),COV)
    d<-mvrnorm(1,rep(0,length(r)),COV2) #mvrnorm(1,rep(0,length(r)),COV_r2)
    N[t+1,] <-N[t,] * exp(r * (1 - (a *N[t,]/KK))+e+d/((N[t,])^0.5))
    N[t+1,][N[t+1,]<=0] <- NaN
  }
  return(N) }

LV_mix <- function(r, a, N0, steps,var_e,var_d) {
  N <- matrix(c(N0, rep(numeric(steps - 1), length(r))), steps, length(r), byrow=T)
  COV <- matrix(rep(-0.1*var_e,sr_i^2),nrow=sr_i,ncol=sr_i)
  diag(COV) <- var_e
  COV2 <- matrix(rep(-0.1*var_d,sr_i^2),nrow=sr_i,ncol=sr_i)
  diag(COV2) <- var_d
  for(t in 1: (steps - 1)) {
    e<-mvrnorm(1,rep(0,length(r)),COV)
    d<-mvrnorm(1,rep(0,length(r)),COV2)#mvrnorm(1,rep(0,length(r)),COV_r2)
    N[t+1,] <-N[t,] * exp(r * (1 - (a %*% N[t,]/KK))+e+d/((N[t,])^0.5))#*rbinom(length(r),1,0.5)
    N[t+1,][N[t+1,]<=0] <- NaN
  }
  return(N) }
##rm1=rm2=1====

var_e<-0.01^2
var_d<-0.2^2
sr <- seq(2, 6, 2)#species richness
sim_com_N <- data.frame()
com <- 1:100 #community每个多样性梯度有100个群落，在一百个群落里面改变生态位差异（通过改变竞争系数???
for (i in 1: length(sr)) {
  for (j in 1:length(com)) {
    
    sr_i <- sr[i]					# number of species
    print(paste("sr", sr_i, "com", j, sp=""))##显示运行到哪??? 物种丰富度sr  群落
    
    
    ## community matrix 竞争系数 生态位差异越大，种间竞争越小。算生态位差异和竞争能力差???
    ##这种方法的竞争是对称的，不改变适合度差异，因为2???1???1???2一样，只看生态位差异???
    mat<- matrix(runif(sr_i^2,0,1),sr_i,sr_i)
    diag(mat)<-1.5
    
    ## calculate niche/fitness differences
    
    n_dif <- matrix(nrow=sr_i, ncol=sr_i)
    f_dif <- matrix(nrow=sr_i, ncol=sr_i)
    
    for (k in 1:sr_i) {
      for (l in k+1:sr_i){
        if(l > sr_i) break()
        
        n_dif[k,l] <- 1- sqrt((mat[k,l]*mat[l,k])/(mat[k,k]*mat[l,l]))
        f_dif[k,l] <- sqrt((mat[k,k]*mat[k,l])/(mat[l,l]*mat[l,k]))
        if(f_dif[k,l]<1) 
          f_dif[k,l]<-1/f_dif[k,l]
      }
    }
    n_dif_m <- mean(n_dif, na.rm=T)##有多个物种求生态位差异，要求平均值。但现在的做法不???,随机产生就不???
    f_dif_m <- mean(f_dif, na.rm=T)
    nf_dif_m <- n_dif_m - f_dif_m
    
    ## calculate structural niche/fitness differences @ method based on Saavedra et al. (2017) EcolMono
    ##结构生态位差异，基于那篇文???
    #str_n_dif <- Omega(mat)
    #str_f_dif <- theta(mat, rep(1, sr_i))
    #str_nf_dif <- str_n_dif - str_f_dif
    
    
    ## simulate biomass at equilibrium
    ## ===========================================
    # set initial biomass for monoculture and mixture
    ####t就是step,稳态的力量，波动由随机因子造成，但波动在平衡值附近波???
    t <- 500 		# time sequence
    r <- rep(1, sr_i) 				# intrinsic growith rate per species
    #r <- runif(sr_i, 0.5, 2) # intrinsic growith rate per species
    ##r不同影响很大，结果会变化
    
    N0_mono <- round(rnorm(sr_i, mean=1, sd=0.1))#取整
    N0_mix <- N0_mono/sr_i#种内竞争大于种间竞争，结果不太受初始值影???
    alpha_mono <- diag(mat)#对角线，种内竞争。单波直接用种内竞争系数
    alpha_mix <- mat#混播要用整个竞争系数矩阵
    
    sim_mono <- LV_mono(r, alpha_mono, N0_mono, t,var_e,var_d) ## simulation for monoculture
    sim_mix <- LV_mix(r, alpha_mix, N0_mix, t,var_e,var_d) ## simulation for mixture
    
    matplot(1:t, sim_mono[, 1:sr_i], type="l", lwd=1.5, xlab="Time", ylab="Biomass", main="Monoculture")
    matplot(1:t, sim_mix[, 1:sr_i], type="l", lwd=1.5, xlab="Time", ylab="Biomass", main="Mixture")	
    ##不同生物多样性的群落种间竞争系数不同，生态位差异不同，竞争系数相???
    bio_mono <- sim_mono[351:t,]##101???300???100之后就达到了稳态，来算稳定???
    bio_mix <- sim_mix[351:t,]
    
    # community stability
    
    com_mean <- mean(rowSums(bio_mix))##群落生物量，再求平均生物???
    com_sd <- sd(rowSums(bio_mix))
    com_stab <- com_mean/com_sd;sum_var<-sum(diag(cov(bio_mix)));sum_cov<-sum(cov(bio_mix))-sum_var
    
    # synchrony
    syn <- data.frame()
    for (n in 1:nrow(bio_mix)) {
      bio_mix_n <- bio_mix[n,]
      sp_i <- paste("sp", 1:sr_i, sep="")
      t_n <- rep(n, sr_i)
      syn_n <- cbind(sp_i, t_n, bio_mix_n)
      syn <- rbind(syn, syn_n)
    }
    
    colnames(syn) <- c("species", "time", "biomass")
    syn$time <- as.numeric(syn$time)
    syn$biomass <- as.numeric(syn$biomass)
    syn_l <- synchrony(syn, time.var="time", species.var="species", abundance.var="biomass")
    #syn_g <- synchrony(syn, time.var="time", species.var="species", abundance.var="biomass", metric="Gross")
    syn_g <- 0
    
    # extinction
    ext <- length(bio_mix[bio_mix<0.0001])##哪些物种消失了，应该看最后一行，改一???
    ext_p <- ext/sr_i#消失的比???
    sr_r <- sr_i - ext
    
    # calculate complementarity/selection Mhichal Lorean
    RYT <- sum(colMeans(bio_mix)/colMeans(bio_mono))#总相对产量，和comp强相???
    dRY <- colMeans(bio_mix)/colMeans(bio_mono)- 1/sr_i
    comp <- sr_i*mean(dRY)*mean(colMeans(bio_mono))#和生态位差异有关???
    sel <- sr_i*cov(dRY, colMeans(bio_mono))#和抽样效应有关系
    delta_Y <- sum(dRY*colMeans(bio_mono)) # or delta_Y=comp+sel  总生物多样性效???
    eveness<-1-2/pi*atan(sum((log(sum(colMeans(bio_mix)))-sum(log(colMeans(bio_mix)))/sr_i)^2)/sr_i)
    syn<-log(var(rowSums(bio_mix))/sum(diag(var(bio_mix))))
    
    comb <- c(sr_i, sum(colMeans(bio_mix)), RYT, delta_Y, comp, sel, com_stab, com_mean, com_sd, syn_l, syn_g, n_dif_m, f_dif_m, nf_dif_m,  ext, ext_p, sr_r,sum_var,sum_cov,eveness,syn)
    sim_com_N <- rbind(sim_com_N, comb)
    ##不同生物多样性的群落本质上遵循同样的规律，物种生态位差异逐渐增长。整体的原则是一百个群落，生态位差异逐渐增长???
  }
}
colnames(sim_com_N) <- c("diversity", "productivity", "RYT", "delta_Y", "comp", "sel", "com_stab", "com_mean", "com_sd", "syn_l", "syn_g", "n_dif_m", "f_dif_m", "nf_dif_m" , "ext", "ext_p", "diversity_r","sum_var","sum_cov","eveness","syn")
spe_sd<-(sim_com_N$sum_var/sim_com_N$diversity)^0.5
spe_mean<-sim_com_N$com_mean/sim_com_N$diversity
spe_stab<-spe_mean/spe_sd
sim_com_N<-cbind(sim_com_N,spe_stab,spe_mean,spe_sd)
colnames(sim_com_N) <- c("diversity", "productivity", "RYT", "delta_Y", "comp", "sel", "com_stab", "com_mean", "com_sd", "syn_l", "syn_g", "n_dif_m", "f_dif_m", "nf_dif_m" , "ext", "ext_p", "diversity_r","sum_var","sum_cov","eveness","syn","spe_stab","spe_mean","spe_sd")
write.table(sim_com_N, "D:/学习/研究???/论文/biodiversity/alphaij/sim_com_N_div_ij.txt", row.names=F)
sim_com_N<-read.table("sim_com_N_div_ij.txt",T)
##plot====
plot(com_stab ~ n_dif_m, col=diversity/2, data=sim_com_N)
legend("topleft", legend=seq(2, 10, 2), col=seq(1, 10, 1), title="Diversity", box.col = "white", pch=1, cex=0.5)
plot(com_sd ~ n_dif_m, col=diversity/2, data=sim_com_N)
legend("topright", legend=seq(2, 10, 2), col=seq(1, 10, 1), title="Diversity", box.col = "white", pch=1, cex=0.5)
plot(com_mean ~ n_dif_m, col=diversity/2, data=sim_com_N)
legend("bottomleft", legend=seq(2, 10, 2), col=seq(1, 10, 1), title="Diversity", box.col = "white", pch=1, cex=0.5)
cor.test(sim_com_N$com_sd, sim_com_N$n_dif_m)

plot(com_stab ~ f_dif_m, col=diversity/2, data=sim_com_N)
legend("topright", legend=seq(2, 10, 2), col=seq(1, 10, 1), title="Diversity", box.col = "white", pch=1, cex=0.5)
plot(com_sd ~ f_dif_m, col=diversity/2, data=sim_com_N)
legend("bottomright", legend=seq(2, 10, 2), col=seq(1, 10, 1), title="Diversity", box.col = "white", pch=1, cex=0.5)
plot(com_mean ~ f_dif_m, col=diversity/2, data=sim_com_N)
legend("bottomright", legend=seq(2, 10, 2), col=seq(1, 10, 1), title="Diversity", box.col = "white", pch=1, cex=0.5)
cor.test(sim_com_N$com_sd, sim_com_N$n_dif_m)

syn<-(sim_com_N$com_sd/(sim_com_N$spe_sd*sim_com_N$diversity))^2
plot(syn ~ n_dif_m, col=diversity/2, data=sim_com_N)
legend("topleft", legend=seq(2, 10, 2), col=seq(1, 10, 1), title="Diversity", box.col = "white", pch=1, cex=0.5)

plot(syn ~ f_dif_m, col=diversity/2, data=sim_com_N,xlim=c(1,5))
legend("topright", legend=seq(2, 10, 2), col=seq(1, 10, 1), title="Diversity", box.col = "white", pch=1, cex=0.5)

cor.test(syn, sim_com_N$n_dif_m)
cor.test(syn, sim_com_N$f_dif_m)


##ρ=±1====
c<-2
phi_e<-0.5
#var<-0.01^2
KK<-20000

LV_mono <- function(r, a, N0, steps,var_e,var_d) {
  N <- matrix(c(N0, rep(numeric(steps - 1), length(r))), steps, length(r), byrow=T)
  for(t in 1: (steps - 1)) { 
    e<-mvrnorm(1,0,var_e)
    d<-mvrnorm(1,0,var_d)#mvrnorm(1,rep(0,length(r)),COV_r2)
    N[t+1,] <-N[t,] * exp(r * (1 - (a *N[t,]/KK))+e*sample(c(1,-1),length(r),T)+d*sample(c(1,-1),length(r),T)/((N[t,])^0.5))
    N[t+1,][N[t+1,]<=0] <- NaN
  }
  return(N) }

LV_mix <- function(r, a,N0, steps,var_e,var_d) {
  N <- matrix(c(N0, rep(numeric(steps - 1), length(r))), steps, length(r), byrow=T)
  for(t in 1: (steps - 1)) {
    e<-mvrnorm(1,0,var_e)
    d<-mvrnorm(1,0,var_d)#mvrnorm(1,rep(0,length(r)),COV_r2)
    N[t+1,] <-N[t,] * exp(r * (1 - (a %*% N[t,]/KK))+e*sample(c(1,-1),length(r),T)+d*sample(c(1,-1),length(r),T)/((N[t,])^0.5))#*rbinom(length(r),1,0.5)
    N[t+1,][N[t+1,]<=0] <- NaN
  }
  return(N) }
##rm1=rm2=1==== ·

var_e<-0.01^2
var_d<-0.2^2
sr <- seq(2,10 , 2)#species richness
sim_com_N <- data.frame()
com <- 1:100 #community每个多样性梯度有100个群落，在一百个群落里面改变生态位差异（通过改变竞争系数???
for (i in 1: length(sr)) {#
  for (j in 1:length(com)) {
    
    sr_i <- sr[i]					# number of species
    print(paste("sr", sr_i, "com", j, sp=""))##显示运行到哪??? 物种丰富度sr  群落
    
    
    ## community matrix 竞争系数 生态位差异越大，种间竞争越小。算生态位差异和竞争能力差???
    ##这种方法的竞争是对称的，不改变适合度差异，因为2???1???1???2一样，只看生态位差异???
    mat<- matrix(runif(sr_i^2,0,1),sr_i,sr_i)
    diag(mat)<-2
    
    ## calculate niche/fitness differences
    
    n_dif <- matrix(nrow=sr_i, ncol=sr_i)
    f_dif <- matrix(nrow=sr_i, ncol=sr_i)
    
    for (k in 1:sr_i) {
      for (l in k+1:sr_i){
        if(l > sr_i) break()
        
        n_dif[k,l] <- 1- sqrt((mat[k,l]*mat[l,k])/(mat[k,k]*mat[l,l]))
        f_dif[k,l] <- sqrt((mat[k,k]*mat[k,l])/(mat[l,l]*mat[l,k]))
        if(f_dif[k,l]<1) 
          f_dif[k,l]<-1/f_dif[k,l]
      }
    }
    n_dif_m <- mean(n_dif, na.rm=T)##有多个物种求生态位差异，要求平均值。但现在的做法不???,随机产生就不???
    f_dif_m <- mean(f_dif, na.rm=T)
    nf_dif_m <- n_dif_m - f_dif_m
    
    ## calculate structural niche/fitness differences @ method based on Saavedra et al. (2017) EcolMono
    ##结构生态位差异，基于那篇文???
    str_n_dif <- Omega(mat)
    str_f_dif <- theta(mat, rep(1, sr_i))
    str_nf_dif <- str_n_dif - str_f_dif
    
    
    ## simulate biomass at equilibrium
    ## ===========================================
    # set initial biomass for monoculture and mixture
    ####t就是step,稳态的力量，波动由随即银子造成，但波动在平衡值附近波???
    t <- 500 		# time sequence
    r <- rep(1, sr_i) 				# intrinsic growith rate per species
    #r <- runif(sr_i, 0.5, 2) # intrinsic growith rate per species
    ##r不同影响很大，结果会变化
    
    N0_mono <- round(rnorm(sr_i, mean=1, sd=0.1))#取整
    N0_mix <- N0_mono/sr_i#种内竞争大于种间竞争，结果不太受初始值影???
    alpha_mono <- diag(mat)#对角线，种内竞争。单波直接用种内竞争系数
    alpha_mix <- mat#混波要用整个竞争系数矩阵
    
    sim_mono <- LV_mono(r, alpha_mono, N0_mono, t,var_e,var_d) ## simulation for monoculture
    sim_mix <- LV_mix(r, alpha_mix, N0_mix, t,var_e,var_d) ## simulation for mixture
    
    matplot(1:t, sim_mono[, 1:sr_i], type="l", lwd=1.5, xlab="Time", ylab="Biomass", main="Monoculture")
    matplot(1:t, sim_mix[, 1:sr_i], type="l", lwd=1.5, xlab="Time", ylab="Biomass", main="Mixture")	
    ##不同生物多样性的群落种间竞争系数不同，生态位差异不同，竞争系数相???
    bio_mono <- sim_mono[351:t,]##101???300???100之后就达到了稳态，来算稳定???
    bio_mix <- sim_mix[351:t,]
    
    # community stability
    
    com_mean <- mean(rowSums(bio_mix))##群落生物量，再求平均生物???
    com_sd <- sd(rowSums(bio_mix))
    com_stab <- com_mean/com_sd;sum_var<-sum(diag(cov(bio_mix)));sum_cov<-sum(cov(bio_mix))-sum_var
    
    # synchrony
    syn <- data.frame()
    for (n in 1:nrow(bio_mix)) {
      bio_mix_n <- bio_mix[n,]
      sp_i <- paste("sp", 1:sr_i, sep="")
      t_n <- rep(n, sr_i)
      syn_n <- cbind(sp_i, t_n, bio_mix_n)
      syn <- rbind(syn, syn_n)
    }
    
    colnames(syn) <- c("species", "time", "biomass")
    syn$time <- as.numeric(syn$time)
    syn$biomass <- as.numeric(syn$biomass)
    syn_l <- synchrony(syn, time.var="time", species.var="species", abundance.var="biomass")
    if(is.na(sim_mix)[500,1]=="TRUE") {syn_g<-0} else {syn_g <- synchrony(syn, time.var="time", species.var="species", abundance.var="biomass", metric="Gross") }
       syn_g<-0 
    
    # extinction
    ext <- length(bio_mix[bio_mix<0.0001])##哪些物种消失了，应该看最后一行，改一???
    ext_p <- ext/sr_i#消失的比???
    sr_r <- sr_i - ext
    
    # calculate complementarity/selection Mhichal Lorean
    RYT <- sum(colMeans(bio_mix)/colMeans(bio_mono))#总相对产量，和comp强相???
    dRY <- colMeans(bio_mix)/colMeans(bio_mono)- 1/sr_i
    comp <- sr_i*mean(dRY)*mean(colMeans(bio_mono))#和生态位差异有关???
    sel <- sr_i*cov(dRY, colMeans(bio_mono))#和抽样效应有关系
    delta_Y <- sum(dRY*colMeans(bio_mono)) # or delta_Y=comp+sel  总生物多样性效???
    eveness<-1-2/pi*atan(sum((log(sum(colMeans(bio_mix)))-sum(log(colMeans(bio_mix)))/sr_i)^2)/sr_i)
    syn<-log(var(rowSums(bio_mix))/sum(diag(var(bio_mix))))
    
    comb <- c(sr_i, sum(colMeans(bio_mix)), RYT, delta_Y, comp, sel, com_stab, com_mean, com_sd, syn_l, syn_g, n_dif_m, f_dif_m, nf_dif_m, ext, ext_p, sr_r,sum_var,sum_cov,eveness,syn)
    sim_com_N <- rbind(sim_com_N, comb)
    ##不同生物多样性的群落本质上遵循同样的规律，物种生态位差异逐渐增长。整体的原则是一百个群落，生态位差异逐渐增长???
  }
}
colnames(sim_com_N) <- c("diversity", "productivity", "RYT", "delta_Y", "comp", "sel", "com_stab", "com_mean", "com_sd", "syn_l", "syn_g", "n_dif_m", "f_dif_m", "nf_dif_m", "ext", "ext_p", "diversity_r","sum_var","sum_cov","eveness","syn")

spe_sd<-(sim_com_N$sum_var/sim_com_N$diversity)^0.5
spe_mean<-sim_com_N$com_mean/sim_com_N$diversity
spe_stab<-spe_mean/spe_sd
sim_com_N<-cbind(sim_com_N,spe_stab,spe_mean,spe_sd)
colnames(sim_com_N) <- c("diversity", "productivity", "RYT", "delta_Y", "comp", "sel", "com_stab", "com_mean", "com_sd", "syn_l", "syn_g", "n_dif_m", "f_dif_m", "nf_dif_m", "ext", "ext_p", "diversity_r","sum_var","sum_cov","eveness","syn","spe_stab","spe_mean","spe_sd")
write.table(sim_com_N, "D:/ѧϰ/????/biodiversity/alphaij/sim_com_N_div_ij.txt", row.names=F)


# Code for reproducing simulation study results published at:

# Vilor-Tejedor N, Gonzalez JR, Calle ML. Efficient and Powerful Method for
# Combining P-Values in Genome-Wide Association Studies. IEEE/ACM Trans Comput Biol
# Bioinform. 2016 Nov;13(6):1100-1106. doi: 10.1109/TCBB.2015.2509977. Epub 2015
# Dec 22. PubMed PMID: 28055892.

# Install and check packages


#library("globalGSA")
library("mc2d")
library("SNPassoc")
library("genetics")
source("/home/nvilor/hapgen2/Code_new.R")


# Keep simulated data files:
files <- dir()

# Function to count the total number of significant p-values:
countSP <- function(pvalors) {
	return(sum(pvalors<=0.05))
}

# Function to generate LDmatrix

LDmatrixGen <- function(data) {
	snps <- setupSNP(data, colSNPs=2:ncol(data), name.genotypes=c(0,1,2))
	snps <- makeGenotypes(data.frame(snps[,-1]))
	ldgsa <- LD(snps)
	rmatrix <- ldgsa$r
	rmatrix[is.na(rmatrix)] <- 0
	rtmatrix <- t(rmatrix)
	LDmatrix <- rmatrix+rtmatrix
	diag(LDmatrix) <- 1
	return(LDmatrix)
}



# Gene-set analysis functions: 
globalARTP_par <- function(i, file, ncores, B, K, output) {
	print(file[i])
	data <- read.table(file[i], header=T)
	core <- (i%%ncores)+1
	a <- proc.time()
	set.seed(123456789)
	ans_global <- globalARTP(data, B, K)
	ans_time <- proc.time()-a
	write.table(ans_global, paste(output,"_output_", core, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE,append=TRUE)
	print(ans_global)
	write.table(ans_time[3], paste(output,"_output_time_", core, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE,append=TRUE)
	return(0)
}


ARTP_par <- function(i, file, ncores, B, K, output) {
	print(file[i])
	data=read.table(files[i], header=T)
	core <- (i%%ncores)+1
	a <- proc.time()
	set.seed(123456789)
	ans_global <- globalARTP(data, B, K, addit=TRUE)
	ans_time <- proc.time()-a
	write.table(ans_global, paste(output,"_output_", core, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE,append=TRUE)
	print(ans_global)
	write.table(ans_time[3], paste(output,"_output_time_", core, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE,append=TRUE)
	return(0)
}

globalEVT_par <- function(i, file, ncores, K, output, LD=FALSE) {
	print(file[i])
	data<-read.table(files[i], header=T)
	core <- (i%%ncores)+1
	if(LD==FALSE) {
		a <- proc.time()
		set.seed(987654321)
		ans_global <- globalEVT(data, K)
		ans_time <- proc.time()-a

	} else {
		LDmat = LDmatrixGen(data)
		a <- proc.time()
		set.seed(987654321)
		ans_global <- globalEVT(data, K, LDmatrix=LDmat)
		ans_time <- proc.time()-a
	}
	write.table(ans_global, paste(output,"_output_", core, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE,append=TRUE)
	print(ans_global)
	write.table(ans_time[3], paste(output,"_output_time_", core, ".txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE,append=TRUE)
	return(0)
}








# Computation:

nCPU <- 45

mclapply(1:length(files), globalEVT_par, file = files, ncores=nCPU, K=5, LD=TRUE, output=paste("GSA_Results_globalEVT_Scenario", substr(files[1],9,9), sep=""), mc.cores=nCPU, mc.preschedule=TRUE)

mclapply(1:length(files), globalARTP_par, file = files, ncores=nCPU, B=1000, K=5, output=paste("GSA_Results_globalARTP_Scenario", substr(files[1],9,9), sep=""), mc.cores=nCPU, mc.preschedule=TRUE)

mclapply(1:length(files), ARTP_par, file = files, ncores=nCPU, B=1000, K=5, output=paste("GSA_Results_ARTP_Scenario", substr(files[1],9,9), sep=""), mc.cores=nCPU, mc.preschedule=TRUE)



# Summary of results:

#globalARTP

out.globalARTP <- NULL
for(i in 1:nCPU) {
	data <- read.table(paste(paste("GSA_Results_globalARTP_Scenario", substr(files[1],9,9), sep=""), "_output_", i, ".txt", sep=""))
	out.globalARTP <- rbind(out.globalARTP, data)
}
names(out.globalARTP) <- c("NPerm", "K", "Kopt", "genevalue")
(res.gsa <- data.frame(Scenario=substr(files[1],9,9),res_globalARTP=countSP(out.globalARTP)))
write.table(res.gsa, paste("GSA_Results_globalARTP_Scenario", substr(files[1],9,9), ".txt", sep=""), sep="\t", row.names=F)

out.globalARTP_time <- NULL
for(i in 1:nCPU) {
	data <- read.table(paste(paste("GSA_Results_globalARTP_Scenario", substr(files[1],9,9), sep=""), "_output_time_", i, ".txt", sep=""))
	out.globalARTP_time <- rbind(out.globalARTP_time, data)
}
names(out.globalARTP_time) <- c( "time")
(res.gsa_time <- data.frame(Scenario=substr(files[1],9,9),res_globalARTP_time=mean(out.globalARTP_time$time, rm.na=T)))
write.table(res.gsa_time, paste("GSA_Results_globalARTP_Scenario_time", substr(files[1],9,9), ".txt", sep=""), sep="\t", row.names=F)



# ARTP

out.ARTP <- NULL
for(i in 1:nCPU) {
	data <- read.table(paste(paste("GSA_Results_ARTP_Scenario", substr(files[1],9,9), sep=""), "_output_", i, ".txt", sep=""))
	out.ARTP <- rbind(out.ARTP, data)
}
names(out.ARTP) <- c("NPerm", "K", "Kopt", "genevalue")
(res.gsa <- data.frame(Scenario=substr(files[1],9,9),res_ARTP=countSP(out.ARTP)))
write.table(res.gsa, paste("GSA_Results_ARTP_Scenario", substr(files[1],9,9), ".txt", sep=""), sep="\t", row.names=F)

out.ARTP_time <- NULL
for(i in 1:nCPU) {
	data <- read.table(paste(paste("GSA_Results_ARTP_Scenario", substr(files[1],9,9), sep=""), "_output_time_", i, ".txt", sep=""))
	out.ARTP_time <- rbind(out.ARTP_time, data)
}
names(out.ARTP_time) <- c( "time")
(res.gsa_time <- data.frame(Scenario=substr(files[1],9,9),res_ARTP_time=mean(out.ARTP_time$time)))
write.table(res.gsa_time, paste("GSA_Results_ARTP_Scenario_time", substr(files[1],9,9), ".txt", sep=""), sep="\t", row.names=F)



#globalEVT

out.globalEVT <- NULL
for(i in 1:nCPU) {
	data <- read.table(paste(paste("GSA_Results_globalEVT_Scenario", substr(files[1],9,9), sep=""), "_output_", i, ".txt", sep=""))
	out.globalEVT <- rbind(out.globalEVT, data)
}
names(out.globalEVT) <- c("K", "genevalue")
out.globalEVT <- out.globalEVT[!is.na(out.globalEVT[,2]),]

(res.gsa <- data.frame(Scenario=substr(files[1],9,9),res_globalEVT=countSP(out.globalEVT)))
write.table(res.gsa, paste("GSA_Results_globalEVT_Scenario", substr(files[1],9,9), ".txt", sep=""), sep="\t", row.names=F)

out.globalEVT_time <- NULL
for(i in 1:nCPU) {
	data <- read.table(paste(paste("GSA_Results_globalEVT_Scenario", substr(files[1],9,9), sep=""), "_output_time_", i, ".txt", sep=""))
	out.globalEVT_time <- rbind(out.globalEVT_time, data)
}
names(out.globalEVT_time) <- c( "time")

(res.gsa_time <- data.frame(Scenario=substr(files[1],9,9),res_globalEVT_time=mean(out.globalEVT_time$time)))
write.table(res.gsa_time, paste("GSA_Results_globalEVT_Scenario_time", substr(files[1],9,9), ".txt", sep=""), sep="\t", row.names=F)



# Together:
res.gsa <- data.frame(Scenario=substr(files[1],9,9),res_globalEVT=countSP(out.globalEVT$genevalue),res_globalARTP=countSP(out.globalARTP$genevalue),res_ARTP=countSP(out.ARTP$genevalue))
write.table(res.gsa, paste("GSA_Results_Scenario", substr(files[1],9,9), ".txt", sep=""), sep="\t", row.names=F)


# Together:
(res.gsa_time <- data.frame(Scenario=substr(files[1],9,9),res_globalEVT_time=mean(out.globalEVT_time$time), res_globalARTP_time=mean(out.globalARTP_time$time),res_ARTP_time=mean(out.ARTP_time$time)))
write.table(res.gsa_time, paste("GSA_time_Results_Scenario", substr(files[1],9,9), ".txt", sep=""), sep="\t", row.names=F)





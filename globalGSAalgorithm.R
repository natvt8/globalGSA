## Biological functions for codifying genotypes:
Codadd <- function(x) {return(as.numeric(factor(x))-1)}
Coddom <- function(x) {return(ifelse(Codadd(x)>0,1,0))}
Codrec <- function(x) {return(ifelse(Codadd(x)>1,1,0))}
Codcodom <- function(x) {return(as.factor(Codadd(x)))}


## Data codification: 
ff <- function(data, model="add") {
	Geno <- data[,2:ncol(data)]
    	if(model=="add") {
		Geno <- data.frame(sapply(Geno, Codadd))		
	}
	if(model=="codom") {
		Geno <- data.frame(sapply(Geno, Codcodom))
	}
	if(model=="dom") {
    		Geno <- data.frame(sapply(Geno, Coddom))
	}
	if(model=="rec") {
    		Geno <- data.frame(sapply(Geno, Codrec))
	}
	return(data.frame(data[,1],Geno))
}


# Create formula:
CreateFormula <- function(data, i, covariable=NULL) {
	if(length(covariable)==0) { fmla <- as.formula(data[,1] ~ data[,i]); return(fmla)
	} else {
		covariable <- as.data.frame(covariable)
		covname <- colnames(covariable)
		fmla <- as.formula(paste("data[,1] ~ data[,i] + ", paste(covname, collapse= "+")))
	}	
	return(fmla)
}


## Generation pvalues:
pvalFmla <- function(data, i, covariable=NULL, family=binomial) {
		if(length(unique(data[,i])) == 1) { 
			pval <- 1
      	} else {
			form=CreateFormula(data,i,covariable)
			if(length(covariable)==0) { 
				fit <- glm(formula=form, data=data, family=family)
				pval <- anova(fit, test="LRT")[2,5]
				return(pval)
			}
			fit <- glm(formula=form, data=data.frame(data, covariable), family=family)
		  	pval <- anova(fit, test="LRT")[2,5]
		}
		return(pval)
}



# Computation of inheritance min pvalues:
runPvalues <- function(data, addit = FALSE, covariable=NULL, family="binomial") {
	padd <- NULL  
	dataA <- ff(data, model="add")
	if(addit==FALSE) {
		dataD <- ff(data, model="dom")
		dataR <- ff(data, model="rec")
		dataC <- ff(data, model="codom")
		pdom <- prec <- pcodom <- NULL
		for(i in 2:ncol(data)){
			pdom <- c(pdom, pvalFmla(dataD, i, covariable, family))
			prec <- c(prec, pvalFmla(dataR, i, covariable, family))
			pcodom <- c(pcodom, pvalFmla(dataC, i, covariable, family))
			padd <- c(padd, pvalFmla(dataA, i, covariable, family))
		}
		pvalors <- data.frame(pdom, pcodom, prec, padd)
    	pvalors$min <- apply(pvalors,1,min,na.rm=TRUE)
    	return(pvalors$min)
	}
	if (addit==TRUE) { 
		for(i in 2:ncol(data)){
			padd <- c(padd, pvalFmla(dataA, i, covariable, family)) 
		}
	} 	
  	return(padd)
}


# Run permutations:  
runPermut <- function(data, addit=FALSE, covariable=NULL, family="binomial") {
	TraitR <- sample(data[,1])
   	dataR <- data.frame(TraitR, data[,2:ncol(data)])
    	return(runPvalues(dataR, addit, covariable, family))
  }


# Generate minPvalues for each permutation:
GeneratePvalues <- function(data, B, addit=FALSE, covariable=NULL, family="binomial") {
  pvalors <- numeric(ncol(data)-1)
  pvalors <- sort(runPvalues(data, addit, covariable,family))
  pvalors <- rbind(pvalors, t(sapply(1:B, function(i) sort(runPermut(data, addit, covariable, family)))))
  rownames(pvalors) <- NULL  
  return(pvalors)
}

# Select genes from an specific dataset:
Selected_genes <- function(Gene, gene_list, data) {
  idX <- as.character(gene_list$Id[gene_list$Gene == Gene])
  return(data[,idX])
}


################ Fisher Combination method ###################

fisher <- function(x) {1-pchisq(-2 * sum(log(x)),df=2*length(x))}

globalFisher <-
  function(data, B, gene_list, Gene="all", addit=FALSE, covariable=NULL, family="binomial") {
    Geno <- data[,-1]
    Trait <- factor(data[,1])
    output <- NULL
    output$nPerm <- B
    if(Gene!="all") {
      data <- data.frame(Trait, Selected_genes(Gene, gene_list, data))
      output$Gene <- Gene  
    }
    x <- GeneratePvalues(data, B, addit, covariable, family)
    Wb <- apply(x,1,fisher)
    output$genevalue <- sum(Wb[1] <= rank(Wb))/(B+1)
    return(output)
  }


################ Simes combination method ###################
simes <- function(x) {min(length(x)*sort(x)/seq(1,length(x), by=1))}

globalSimes <-
  function(data, B, gene_list, Gene="all", addit=FALSE, covariable=NULL, family=binomial) {
    Geno <- data[,-1]
    Trait <- factor(data[,1])
    output <- NULL
    output$nPerm <- B
    if(Gene!="all") {
      data <- data.frame(Trait, Selected_genes(Gene, gene_list, data))
      output$Gene <- Gene  
    }
    x <- GeneratePvalues(data, B, addit, covariable, family)
    Wb <- apply(x,1,simes)
    output$genevalue <- sum(Wb[1] <= rank(Wb))/(B+1)	
    return(output)
  }

####################### globalARTP ############################

# Truncation process:    
Trunkpoint <- function(B, K, pvalues) {
  trunks <- matrix(0,ncol=K, nrow=B+1)
  trunks[,1] <- -log(pvalues[,1])
  if(K==1) return(trunks)
  for(j in 2:K) {
    trunks[,j] <- trunks[,j-1]+(-log(pvalues[,j]))
  }			
  return(trunks)
}

#Adaptive process:
EstimatePvalue <- function(B, K, pvalues) {
  WbK <- Trunkpoint(B, K, pvalues)
  sb <- matrix(0,ncol=K, nrow=B+1)
  for(j in 1:K) {
    for(i in 1:(B+1)) {
      sb[i,j] <- sum(WbK[,j] >= WbK[i,j])/(B+1)
    }
  }  
  return(sb)
}	

# General function:
globalARTP <- 
 function(data, B, K, gene_list, Gene="all", addit=FALSE,  covariable=NULL, family="binomial") {
  output <- NULL
  output$nPerm <- B
  if(Gene!="all") {
    colnames(gene_list) <- c("Id", "Gene")
    data <- data.frame(data$y, Selected_genes(Gene, gene_list, data))
    output$Gene <- Gene
  }
  Geno <- data[,-1]
  Trait <- factor(data[,1])
  if(K==1 || ncol(Geno) == 1) {
    output$Trunkpoint <- 1
    output$Kopt <- 1
    output$genevalue <- runPvalues(data, addit, covariable, family) 	
    return(output)
  }else{
    if(K>ncol(Geno)) K <- ncol(Geno)
    output$Trunkpoint <- K
    pvalors <- GeneratePvalues(data, B, addit, covariable, family)
    sb <- EstimatePvalue(B,K,pvalors)
    mP <- apply(sb,1,min)
    output$Kopt <- order(sb[1,])[1]
    output$genevalue <- sum(mP <= mP[1])/(B+1)
  }
  return(output)
}



####################### globalEVT ############################

# Function to calculate TrunkStats:
trunkStat <- function(K, pvalues) {
  trunks <- NULL
  trunks[1] <- -log(pvalues[1])
  if(K==1) return(trunks*2)
  for(j in 2:K) {
    trunks[j] <- trunks[j-1]+(-log(pvalues[j]))
  }    	
  return(trunks*2)
}


globalFun <- function(dades, K, addit=FALSE, covariable=NULL, family=binomial, LDmatrix=NULL) {
  #Obtain min-pvalues:
  pmins <- p.adjust(runPvalues(dades, addit, covariable, family),"fdr")
  # Transform p-values:
  pu <- pbeta(pmins,1,2.2) 
  #Sort p-values:
  pu <- sort(pu)
  dades <- apply(dades[,-1], 2, as.integer)
 
  # UNCORRELATED case using the total number of tests:
  tj <- NULL
  for(i in 1:length(pu)) {
    	tj[i] <- pbeta(pu[i], order(pu[i]), length(pu)-order(pu[i])+1)
  }
  if(length(LDmatrix)==0) { return(trunkStat(K, tj))
  } else { 
  # Compute Me (effective number of tests)
	lambdas <- eigen(LDmatrix)$values 
	M <- ncol(dades[,-1])
	Meff <- M-sum(as.numeric(lambdas>1)*(lambdas-1))
	return(trunkStat(Meff,tj))
  } 
}

LDextra <- function(n) {
  if(n <= 1 & n >= 0) { return(n*(3.25+0.75*n))
  } else { return(n*(3.27+0.71*n)) 
  }
}

adaptation <- function(data, K, addit=FALSE, covariable=NULL, family=binomial, LDmatrix=NULL) {
  Uk <- NULL
  if((ncol(data)-1) < K) K = ncol(data)-1
  Statk <- globalFun(data, K=K, addit, covariable, family, LDmatrix)
  if(length(LDmatrix)!=0) {
    LDmatrix <- LDmatrix[lower.tri(LDmatrix, diag = FALSE)==TRUE]
    f <- 2*(2*(ncol(data)-1))^2/(4*(ncol(data)-1)+2*sum(sapply(data,LDextra)))
    c<- 4*(ncol(data)-1)+2*sum(sapply(LDmatrix,LDextra))/2*(2*(ncol(data)-1))
    Uk <- sapply(1:K, function(x, st = Statk/c) 1-pchisq(st, df = f))
    return(min(Uk))
  } else {
    Uk <- sapply(1:K, function(x, st = Statk) 1-pchisq(st, df = 2*x))
    return(min(Uk))    
  } 
}

 

globalEVT <- function(data, K, gene_list, Gene="all", addit=FALSE, covariable=NULL, family=binomial, LDmatrix=NULL) {
	output <- NULL
  	if(Gene != "all") {
		colnames(gene_list) <- c("Id", "Gene")
      	data <- data.frame(data$y, Selected_genes(Gene, gene_list, data))
      	output$Gene <- Gene
  	}
  	if(K==1 || ncol(data) <= 2) {
  		output$Trunkpoint <- 1
    		output$genevalue <- runPvalues(data, addit, covariable, family) 	
    		return(output)
  	}else{
		if(K >(ncol(data)-1)) K <- ncol(data)
		output$Trunkpoint <- K
    		W <- NULL
		dades <- data
		W[1] <- adaptation(data, K)
		for(i in 2:100) {
			dades[,1] <- sample(data[,1])
			W[i] <- adaptation(dades, K, LDmatrix=LDmatrix)
		}
		zeff <- (1-mean(W))/mean(W)
		output$genepvalue <- pbeta(W[1], 1, zeff)
		return(output)
 	}
}

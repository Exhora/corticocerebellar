require(CCA)

load("../dados/phenotypeComplete.RData")                                        
load("../dados/roisInteresse.RData")                                            
load("../dados/allRoisCerebelocc400.RData")
load("../dados/volremovidos.RData")
load("../dados/funcIdade.RData")
     
offset <-2
id <- phenotype$SUB_ID

grupoTeste <- which(funcIdade < 31 & vol.removidos < 5)
corRois <- matrix(0, length(grupoTeste), length(roisCerebro))   
pessoas <- list()
p<-0
roisCerebeloExcluidas <- matrix(0, length(grupoTeste))
roisCerebroExcluidas <- matrix(0, length(grupoTeste))
for(ind in grupoTeste){
	#carrega o arquivo da pessoa id[ind]
	file<-paste0(paste0("../original/pm.ORIGINAL.00",id[ind], sep=""),".cleanEPI_cc400_TCs.1D", sep="")
	#file<-paste0(paste0("../original/pm.ORIGINAL.00",id[ind], sep=""),".cleanEPI_aal_TCs.1D", sep="")
	if(file.exists(file)){
		print(ind)
	   	p<-p+1
		pessoas[[p]] <- ind
		dadosOrig <- read.table(file, header=T)
		nSeries <- nrow(dadosOrig)
		nRois <- ncol(dadosOrig)-offset
		#pega as series de todas as rois
		dados = matrix(unlist(dadosOrig[ , (offset + 1):(nRois + offset)]), nSeries, nRois, byrow = FALSE)
		
		roisCerebelo <- which(allRoisCerebelo == roisInteresse & roisInteresse == TRUE)
		roisCerebro <- which(allRoisCerebelo != roisInteresse & roisInteresse == TRUE)
		
		seriesCerebelo <- matrix(0, nSeries, length(roisCerebelo))
		seriesCerebro <- matrix(0, nSeries, length(roisCerebro))
		
		#separa em rois do cerebro e rois do cerebelo
		for (i in 1:nSeries){                                                   
			seriesCerebelo[i,] <- dados[i, roisCerebelo]
			seriesCerebro[i,] <- dados[i, roisCerebro]
		}  
		#elimina as ROIs nulas no cerebelo
		seriesCerebelo <- seriesCerebelo[,which(colSums(abs(seriesCerebelo)) > 1)]
		roisCerebeloExcluidas[ind] <- length(which(colSums(abs(seriesCerebelo)) < 1))

		#normaliza as ROIs para variancia 1
		for(i in 1:ncol(seriesCerebelo)){
			dp <- sd(seriesCerebelo[,i])
			seriesCerebelo[,i] <- seriesCerebelo[,i]/dp
		}
		colinear <- list()
		xcor <- cor(seriesCerebelo)
		k<-1
		#verifica colinearidade no cerebelo
		for(i in  2:(ncol(seriesCerebelo)-1)){
			for(j in (i+1):(ncol(seriesCerebelo)-1)){
				if(xcor[i,j] > 0.8){
					colinear[[k]] <- i
					k<-k+1
				}
			}
		}
		#Remover as ROIs colineares
		index <- 1:ncol(seriesCerebelo)
		ROIs <- subset(index, !(index %in% colinear))
		seriesCerebelo <- seriesCerebelo[, ROIs]

		contRoiFora <- 0
		#passa por cada ROI do cerebro, separa as series dessa ROI
		#e faz a correlacao da ROI com todas as ROIs do cerebelo
		serieRoi <- matrix(0, nrow(seriesCerebro), 1)
		for(i in 1:length(roisCerebro)){
			#so faz nas ROIs nao nulas
			if(sum(abs(seriesCerebro[,i]))>1){
				dp <- sd(seriesCerebro[,i])                                     
				serieRoi[,1] <- seriesCerebro[,i]/dp  
				corRois[p,i] <-cc(serieRoi, seriesCerebelo)$cor
				#corRois[p,i] <-gGranger(serieRoi, seriesCerebelo, 1)$B
			}else{
				corRois[p,i] <- NA
				contRoiFora<-contRoiFora+1
			}
		}
		roisCerebroExcluidas[ind] <- contRoiFora
	}
}

pval <- matrix(0, length(roisCerebro))
for(i in 1:length(roisCerebro)){
	pval[i] <- summary(lm(corRois[,i] ~ phenotype$DX_GROUP[grupoTeste]))$coeff[8]
	#pval[i] <- summary(lm(corRois[,i] ~ phenotype$DX_GROUP[grupoTeste] + funcIdade[grupoTeste] + phenotype$SEX[grupoTeste] + factor(phenotype$SITE_ID[grupoTeste])))$coeff[68]
}

adjustp <- p.adjust(pval, method = "fdr", length(pval))
which(adjustp < 0.05)


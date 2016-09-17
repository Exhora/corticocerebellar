source("gGranger.R")

require(CCA)

load("../dados/phenotypeComplete.RData")                                        
load("../dados/roisInteresseAAL.RData")
load("../dados/roisCerebeloAAL.RData")
load("../dados/volremovidos.RData")
load("../dados/funcIdade.RData")
     
offset <-2
id <- phenotype$SUB_ID

roisInteresse <- !roisInteresse
allRoisCerebelo <- roisCerebelo
roisCerebelo <- which(allRoisCerebelo == roisInteresse & roisInteresse == TRUE)
roisCerebro <- which(allRoisCerebelo != roisInteresse & roisInteresse == TRUE)


#grupoTodos <- which(vol.removidos < 5)
grupoTeste <- which(funcIdade < 31 & vol.removidos < 5)
#grupoPaper <- which(funcIdade > 8 & funcIdade < 18 & vol.removidos < 5)
#grupoAdulto <- which(funcIdade > 20 & vol.removidos < 5)
#grupoCriança <- which(funcIdade < 14 & vol.removidos < 5)
#grupoAdolescente <- which(funcIdade > 12 & funcIdade < 17 & vol.removidos < 5)
#grupoJovem <- which(funcIdade > 17 & funcIdade < 25 & vol.removidos < 5)

#grupoTeste <- grupoPaper

corRois <- matrix(0, length(grupoTeste), length(roisCerebro))   
pessoas <- list()
p<-0
ROIsCerebelo <- matrix(0, length(grupoTeste))
roisCerebeloExcluidas <- matrix(0, length(grupoTeste))
roisCerebroExcluidas <- matrix(0, length(grupoTeste))

for(ind in grupoTeste){
	#carrega o arquivo da pessoa id[ind]
	file<-paste0(paste0("../original/pm.ORIGINAL.00",id[ind], sep=""),".cleanEPI_aal_TCs.1D", sep="")
	if(file.exists(file)){
		print(ind)
	   	p<-p+1
		pessoas[[p]] <- ind
		dadosOrig <- read.table(file, header=T)
		nSeries <- nrow(dadosOrig)
		nRois <- ncol(dadosOrig)-offset
		#pega as series de todas as rois
		dados = matrix(unlist(dadosOrig[ , (offset + 1):(nRois + offset)]), nSeries, nRois, byrow = FALSE)
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
			for(j in (i+1):(ncol(seriesCerebelo))){
				if(xcor[i,j] > 0.8){
					if(sum(xcor[i,]) > sum(xcor[j,])){
						colinear[[k]] <- i
					}else{
						colinear[[k]] <- j
					}
					k<-k+1
				}
			}
		}
		
		#Remover as ROIs colineares
		index <- 1:ncol(seriesCerebelo)
		ROIs <- subset(index, !(index %in% colinear))
		seriesCerebelo <- seriesCerebelo[, ROIs]
		ROIsCerebelo[ind] <- length(ROIs) 
		
		contRoiFora <- 0
		#passa por cada ROI do cerebro, separa as series dessa ROI
		#e faz a correlacao da ROI com todas as ROIs do cerebelo
		serieRoi <- matrix(0, nrow(seriesCerebro), 1)
		for(i in 1:length(roisCerebro)){
			#so faz nas ROIs nao nulas
			if(sum(abs(seriesCerebro[,i]))>1){
				dp <- sd(seriesCerebro[,i])                                     
				serieRoi[,1] <- seriesCerebro[,i]/dp  
				#corRois[p,i] <-cc(serieRoi, seriesCerebelo)$cor
				corRois[p,i] <- gGranger(serieRoi, seriesCerebelo, 1)$B
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
	#pval[i] <- summary(lm(corRois[,i] ~ phenotype$DX_GROUP[grupoTeste]))$coeff[8]
	pval[i] <- summary(lm(corRois[,i] ~ factor(phenotype$DX_GROUP[grupoTeste]) + funcIdade[grupoTeste] + phenotype$SEX[grupoTeste] + factor(phenotype$SITE_ID[grupoTeste])))$coeff[68]
	#summary(lm(corRois[,1] ~ phenotype$DX_GROUP[grupoTeste] + funcIdade[grupoTeste] + phenotype$SEX[grupoTeste]))$coeff[14]
}

adjustp <- p.adjust(pval, method = "fdr", length(pval))
signif <- which(adjustp < 0.05)
signif


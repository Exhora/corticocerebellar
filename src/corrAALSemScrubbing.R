#Verificar se a leitura do cerebelo esta completa

require(MASS)

load("../dados/phenotypeComplete.RData")                                        
load("../dados/roisInteresseAAL.RData")
load("../dados/roisCerebeloAAL.RData")
load("../dados/volremovidos.RData")
load("../dados/funcIdade.RData")
     
offset <-2
id <- phenotype$SUB_ID

#removendo o ventriculo
roisInteresse <- !roisInteresse
allRoisCerebelo <- roisCerebelo
roisCerebelo <- which(allRoisCerebelo == roisInteresse & roisInteresse == TRUE)
roisCerebro <- which(allRoisCerebelo != roisInteresse & roisInteresse == TRUE)

#pega so ate 31 anos e menos de 5% de vol removidos
grupoTeste <- which(funcIdade < 31 & vol.removidos < 5)


corRois <- matrix(0, length(grupoTeste), length(roisCerebro))   
corRois2 <- matrix(0, length(grupoTeste), length(roisCerebro))   
p<-0
pessoas <- matrix(0, length(grupoTeste))
acumulada <- matrix(0, length(grupoTeste))
pcs <- matrix(0, length(grupoTeste))
coefcerebelo <- array(0, dim=c(length(grupoTeste), length(roisCerebro), length(roisCerebelo)))

for(ind in grupoTeste){
	#carrega o arquivo da pessoa
	file<-paste0(paste0("../original_sempower/sfnwmrda00",phenotype$SUB_ID[ind], sep=""),"_session_1_rest_1_aal_TCs.1D", sep="")
	print(ind)

	#pega as series de todas as rois
	dadosOrig <- read.table(file, header=T)
	nSeries <- nrow(dadosOrig)
	nRois <- ncol(dadosOrig)-offset
	dados = matrix(unlist(dadosOrig[ , (offset + 1):(nRois + offset)]), nSeries, nRois, byrow = FALSE)
	seriesCerebelo <- matrix(0, nSeries, length(roisCerebelo))      
	seriesCerebro <- matrix(0, nSeries, length(roisCerebro))

	#separa em rois do cerebro e rois do cerebelo
	for (i in 1:nSeries){                                                   
		seriesCerebelo[i,] <- dados[i, roisCerebelo]
		seriesCerebro[i,] <- dados[i, roisCerebro]
	}  
	#normaliza as ROIs para media 0 variancia 1
	for(i in 1:ncol(seriesCerebelo)){
		dp <- sd(seriesCerebelo[,i])
		seriesCerebelo[,i] <- seriesCerebelo[,i]/dp
	}
	#remove as pessoas que tem ROIs zeradas no cerebelo
	if(length(which(is.na(seriesCerebelo)==TRUE))>0) next
	#guarda quantas e quais pessoas efetivamente entraram na analise
	p<-p+1
	pessoas[p] <- ind

		
	#PCA
	res <- prcomp(seriesCerebelo)

	#pega quantos PCs somam 95% para a pessoa p
	pcs[p] <- max(which(cumsum(res$sdev^2/sum(res$sdev^2)) < 0.95)) +1

	serieRoi <- matrix(0, nrow(seriesCerebro), 1)   
	#Faz a regressao das ROIs do cortex com os PCs do cerebelo
	for( i in 1:length(roisCerebro)){
		#verifica se a ROI nao esta zerada
		if(sum(abs(seriesCerebro[,i]))>1){
			#normaliza para media 0, var 1
			dp <- sd(seriesCerebro[,i])                                     
			serieRoi[,1] <- seriesCerebro[,i]/dp
			
			#Faz a regressao com o numero de PCs que somam 95%
			if(pcs[p] == 4)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4]))
			else if(pcs[p] == 5)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5]))
			else if(pcs[p] == 6)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]))
			else if(pcs[p] == 7)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]))
			else if(pcs[p] == 8)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]))
			else if(pcs[p] == 9)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9]))
			else if(pcs[p] == 10)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10]))
			else if(pcs[p]  == 11)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]))
			else if(pcs[p]  == 12)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]))
			else if(pcs[p]  == 13)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]+ res$x[,13]))
			else if(pcs[p]  == 14)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]+ res$x[,13]+ res$x[,14]))
			else if(pcs[p]  == 15)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]+ res$x[,13]+ res$x[,14] + res$x[,15]))
			else if(pcs[p]  == 16)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]+ res$x[,13]+ res$x[,14] + res$x[,15] + res$x[,16]))
			else if(pcs[p]  == 17)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]+ res$x[,13]+ res$x[,14] + res$x[,15] + res$x[,16] + res$x[,17]))
			else if(pcs[p]  == 18)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]+ res$x[,13]+ res$x[,14] + res$x[,15] + res$x[,16] + res$x[,17] + res$x[,18]))
			else if(pcs[p]  == 19)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]+ res$x[,13]+ res$x[,14] + res$x[,15] + res$x[,16] + res$x[,17] + res$x[,18] + res$x[,19]))
			else if(pcs[p]  == 20)
				sm <- summary(lm(serieRoi ~ res$x[,1] + res$x[,2] + res$x[,3] + res$x[,4] + res$x[,5] + res$x[,6]+ res$x[,7]+ res$x[,8]+ res$x[,9] + res$x[,10] + res$x[,11]+ res$x[,12]+ res$x[,13]+ res$x[,14] + res$x[,15] + res$x[,16] + res$x[,17] + res$x[,18] + res$x[,19] + res$x[,20]))

			#salva as correlacoes
			corRois[p,i] <- sm$adj.r.squared

			#salva os coeficientes do cerebelo: matriz de rotacao * coeficientes dos PCs
			coefcerebelo[p,i, ] <- res$rotation[, 1:pcs[p]] %*% as.matrix(sm$coeff[2:(pcs[p]+1)],pcs[p],1)

		}else{
			#se a ROI estiver zerada, ele salva como NA para ser ignorado na regressao
			corRois[p,i] <- NA
		}
	
	}
}

#grupo eh efetivamente entrou na regressao
grupo <- pessoas[1:p]

#coloca o grupo DT como referencia
diagnotico <- phenotype$DX_GROUP[grupo]
diagnotico[which(diagnotico ==2)] = 0

pval_ss <- matrix(0, length(roisCerebro))
for(i in 1:length(roisCerebro)){
	pval_ss[i] <- summary(lm(corRois[1:p,i] ~ diagnotico + funcIdade[grupo] + phenotype$SEX[grupo] + factor(phenotype$SITE_ID[grupo])))$coeff[68]
	#pval_ss[i] <- summary(lm(corRois[1:p,i] ~ factor(phenotype$DX_GROUP[grupo]) + funcIdade[grupo] + phenotype$SEX[grupo] + factor(phenotype$SITE_ID[grupo]) + vol.removidos[grupo]))$coeff[71]
}

signifss <- which(p.adjust(pval_ss, method = "bonferroni", length(pval_ss))<0.05)

#ver o resultado das regressões
sm <- list()
j <- 1
for(i in signif){
	sm[[j]] <- summary(lm(corRois[,i] ~ factor(phenotype$DX_GROUP[grupoTeste]) + funcIdade[grupoTeste] + phenotype$SEX[grupoTeste] + factor(phenotype$SITE_ID[grupoTeste])))
	j<-j+1
}

#ver o beta das regressões (ta errado)
beta <- matrix(0, length(roisCerebro))
for (i in signif){
	beta[i] <- summary(lm(corRois[,i] ~ factor(phenotype$DX_GROUP[grupoTeste]) + funcIdade[grupoTeste] + phenotype$SEX[grupoTeste] + factor(phenotype$SITE_ID[grupoTeste])))$coeff[2]
	print(beta[i])
}


#pegar os betas negativos (mais conectividade no autismo)
betaNegativo <- which(beta<0)
betapositivo <- which(beta>0)

#pegar os betas positivos (menos conectividade no autismo)
roiPositiva <- signif[betapositivo]
roiNegativa <- signif[betaNegativo]

#calcular a media das correlações
media <- matrix(0,length(roisCerebro))
for(i in 1:length(roisCerebro)){
 media[i] <- mean(corRois[,i], na.rm=TRUE)
}

#teste wilcox
wilcoxtest <- matrix(0, length(signif))
j <- 1
for(i in signif){
	wilcoxtest[j] <- wilcox.test(corRois[,i] ~ factor(phenotype$DX_GROUP[grupoTeste]))[3]
	j<-j+1
}

#separa os grupos
autismo <- which(phenotype$DX_GROUP[grupoTeste] == 1)
controle <- which(phenotype$DX_GROUP[grupoTeste] == 2)

#boxplot das regressões
grupo <- phenotype$DX_GROUP[grupoTeste]
grupo[which(grupo == 1)] <- "TEA"
grupo[which(grupo == 2)] <- "DT"

par(mfrow=c(2,3))
boxplot(corRois[,signif[1]] ~ grupo, main = "Right Precentral Gyrus", ylab = "ROI and Cerebellum Correlation", xlab=bquote("Wilcox Test p-value: " ~.(wilcoxtest[1])))
boxplot(corRois[,signif[2]] ~ grupo, main = "Right Precentral Gyrus", ylab = "ROI and Cerebellum Correlation", xlab=bquote("Wilcox Test p-value:"))
boxplot(corRois[,signif[3]] ~ grupo, main = "Right Precentral Gyrus", ylab = "ROI and Cerebellum Correlation", xlab=bquote("Wilcox Test p-value:"))
boxplot(corRois[,signif[4]] ~ grupo, main = "Right Precentral Gyrus", ylab = "ROI and Cerebellum Correlation", xlab=bquote("Wilcox Test p-value:"))
boxplot(corRois[,signif[5]] ~ grupo, main = "Right Precentral Gyrus", ylab = "ROI and Cerebellum Correlation", xlab=bquote("Wilcox Test p-value:"))
boxplot(corRois[,signif[6]] ~ grupo, main = "Right Precentral Gyrus", ylab = "ROI and Cerebellum Correlation", xlab=bquote("Wilcox Test p-value:"))

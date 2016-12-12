source("corrAALSemScrubbing.R")
source("corrAALComScrubbing.R")

quality_pass <- which(abs(pval - pval_ss) < 0.05)

adjustp <- p.adjust(pval[quality_pass], method = "bonferroni", length(quality_pass))
signifp <- adjustp[which(adjustp < 0.05)]
signifRois <- roisCerebro[quality_pass[which(adjustp < 0.05)]]
signif <- quality_pass[which(adjustp < 0.05)]

sm <- list()
j <- 1
for(i in signif){
	sm[[j]] <- summary(lm(corRois[,i] ~ factor(phenotype$DX_GROUP[grupoTeste]) + funcIdade[grupoTeste] + phenotype$SEX[grupoTeste] + factor(phenotype$SITE_ID[grupoTeste])))
	j<-j+1
}


#faz a media dos coeficientes do cerebelo de toda a pop
coefsPop <- matrix(0, length(signif), 26)
coefsAutismo <- matrix(0, length(signif), 26)
coefsControle <- matrix(0, length(signif), 26)
autismo <- which(phenotype$DX_GROUP[pessoas[1:p]] == 1)
controle <- which(phenotype$DX_GROUP[pessoas[1:p]] == 2)
s <- 1
for(j in signif){
	for(i in 1:26){
		coefsPop[s, i] <- mean(coefcerebelo[1:p,j,i])
		coefsAutismo[s,i] <- mean(coefcerebelo[autismo,j,i])
		coefsControle[s,i] <- mean(coefcerebelo[controle,j,i])
	}

	s <- s+1
}

wilcoxon <- matrix(0, 26, 1)

for(i in 1:26)  wilcoxon[i] <- wilcox.test(coefcerebelo[autismo, 1, i], coefcerebelo[controle, 1, i])$p.value
diffcerebelo <- qnorm(1 - wilcoxon, 0, 1)
save(diffcerebelo, file = "diffcerebelo.RData")



                                                                                
rm(list = ls())                                                                 
                                                                                
require(AnalyzeFMRI)                                                            
                                                                                
#carregar aqui o arquivo que tem as rois da regiao a selecionar                 
load('Rois/roisCerebelo.RData')
load('../src/coefsAutismo.RData')
load('../src/coefsControle.RData')

select <- roisCerebelo                                                                                                            
# abrir um arquivo exemplo do dataset, para ligar corretamente cada ROI com seu número original no atlas.
dados = read.table('../original/pm.ORIGINAL.0051212.cleanEPI_aal_TCs.1D', header = T)       
                                                                                
nomes = names(dados)[3:ncol(dados)]                                             
for(i in 1:length(nomes)){                                                      
	  nomes[i] = strsplit(nomes[i], '_')[[1]][2]                                    
}                                                                               
nomes = as.numeric(nomes)                                                       
                                                                                
nomes <- nomes[select]                                                          
                                                                                
#ATLAS = f.read.volume('ADHD200_parcellate_400a.nii')                           
ATLAS = f.read.volume('../templates/aal_mask_pad_Alternativo.nii')                             
                                                                                
SAIDA = array(-2000, dim(ATLAS))
for(i in 1:length(nomes)) SAIDA[which(ATLAS == nomes[i])] = coefsControle[1,i]*10000

h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')                          
f.write.nifti(SAIDA, paste('', paste('controleRoi2', sep = '_'), sep = ''), size = "float", h, nii=TRUE)

SAIDA = array(-2000, dim(ATLAS))
for(i in 1:length(nomes)) SAIDA[which(ATLAS == nomes[i])] = coefsControle[2,i]*10000

h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')                          
f.write.nifti(SAIDA, paste('', paste('controleRoi41', sep = '_'), sep = ''), size = "float", h, nii=TRUE)

SAIDA = array(-2000, dim(ATLAS))
for(i in 1:length(nomes)) SAIDA[which(ATLAS == nomes[i])] = coefsControle[3,i]*10000

h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')                          
f.write.nifti(SAIDA, paste('', paste('controleRoi59', sep = '_'), sep = ''), size = "float", h, nii=TRUE)

SAIDA = array(-2000, dim(ATLAS))
for(i in 1:length(nomes)) SAIDA[which(ATLAS == nomes[i])] = coefsAutismo[1,i]*10000

h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')                          
f.write.nifti(SAIDA, paste('', paste('autismoRoi2', sep = '_'), sep = ''), size = "float", h, nii=TRUE)
SAIDA = array(-2000, dim(ATLAS))
for(i in 1:length(nomes)) SAIDA[which(ATLAS == nomes[i])] = coefsAutismo[2,i]*10000

h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')                          
f.write.nifti(SAIDA, paste('', paste('autismoRoi41', sep = '_'), sep = ''), size = "float", h, nii=TRUE)

SAIDA = array(-2000, dim(ATLAS))
for(i in 1:length(nomes)) SAIDA[which(ATLAS == nomes[i])] = coefsAutismo[3,i]*10000

h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')                          
f.write.nifti(SAIDA, paste('', paste('autismoRoi59', sep = '_'), sep = ''), size = "float", h, nii=TRUE)


load("../src/zvalue.RData")
SAIDA = array(-2000, dim(ATLAS))
for(i in 1:length(nomes)) SAIDA[which(ATLAS == nomes[i])] = zvalue[i]*100

h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')                          
f.write.nifti(SAIDA, paste('', paste('zvalue', sep = '_'), sep = ''), size = "float", h, nii=TRUE)

load("../src/diffcerebelo.RData")
SAIDA = array(-3000, dim(ATLAS))                                                
for(i in 1:length(nomes)) SAIDA[which(ATLAS == nomes[i])] = diffcerebelo[i]*100       
                                                                                
h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')            
f.write.nifti(SAIDA, paste('', paste('diffcerebelo85', sep = '_'), sep = ''), size = "float", h, nii=TRUE)


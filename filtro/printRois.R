                                                                                
rm(list = ls())                                                                 
                                                                                
require(AnalyzeFMRI)                                                            
                                                                                
#carregar aqui o arquivo que tem as rois da regiao a selecionar                 
load('Rois/roisPrint.RData')                                                   
	select <- roisPrint                                                          
                                                                                
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
                                                                                
SAIDA = array(-200, dim(ATLAS))                                                   
for(i in 1:length(nomes)){                                                      
	                #quais voxels pertencem a qual ROI                                      
	          SAIDA[which(ATLAS == nomes[i])] = i
}                                                                               


# as duas proximas linhas gerarao um arquivo clustMap.nii que deve ser aberto no MRIcron.
# abra o MRIcron, vá no menu "File/Open templates/ch2better.nii.gz" (ou ch2bet.nii.gz)
# após, vá no menu "Overlay/add overlay" e abra a imagem gerada (clustMap.nii) na pasta de destino definida
# escolha um mapa de cores adequado (onde estiver escrito "Grayscale" e pronto).
# Voce tambem pode gerar a imagem 3D do cérebro no menu "Window/render" e fatias em "Window/multislice".
h = f.read.nifti.header('../templates/aal_mask_pad_Alternativo.nii')                          
f.write.nifti(SAIDA, paste('', paste('clustMap', sep = '_'), sep = ''), size = "float", h, nii=TRUE)
#versao 1    

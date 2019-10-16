
library(stringr)
library(RColorBrewer)
library(data.table)
library(scales)
library(ggpubr)
library(ggplot2)

af <- read.table("~/oyster/analysis/filtered_allele_freqs.txt", header=TRUE)

#############################################
#############################################
##
## PCA
##
#############################################
#############################################

pops <- colnames(af)[2:7]


pcaResult <- prcomp(t(af[,2:7]), scale=FALSE)

#sal_sig <- which((pcaResult$rotation[,3] > quantile(pcaResult$rotation[,3], 0.95)))

out <- as.data.frame(matrix(ncol=3, nrow=nrow(af)))
colnames(out) <- c("LL", "OP", "TB")
out$LL <- abs(af$LL_33 - af$LL_05)
out$OP <- abs(af$OP_33 - af$OP_05)
out$TB <- abs(af$TB_33 - af$TB_05)


# get proportion of total variance explained:
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(id=pops, Population=substr(pops, 1,2), 
    Treatment=substr(pops, 4,5),
    PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

a <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Population)) +
        geom_point(size=4, color="black", alpha=0.8) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('brown3','blue3'))+
        #theme(legend.position = c(0.88,0.17))+
       # theme(legend.text=element_text(size=8),legend.title=element_blank())+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
       guides(fill=guide_legend(override.aes=list(shape=c(21,21), 
            size=c(5,5), fill=c('brown3','blue3'))))
       #ggtitle("PC1 + PC2")
        #        shape=FALSE,
        #        size=FALSE)+
        #scale_size_manual(values=c(7,5))
#a


dat.p <- data.frame(id=pops, Population=substr(pops, 1,2), 
    Treatment=substr(pops, 4,5),
    PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,3])

a2 <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Population)) +
        geom_point(size=4, color="black", alpha=0.8) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC3: ",percentVar[3],"% variance")) +
        theme_bw() +
       # ylim(-30, 23) + xlim(-50, 65)+
        scale_shape_manual(values=c( 21,22,24))+
        #scale_color_manual(values=c('brown3','blue3', 'darkgoldenrod3', "darkorchid3"))+
        scale_fill_manual(values=c('brown3','blue3'))+
        #theme(legend.position = c(0.88,0.17))+
       # theme(legend.text=element_text(size=8),legend.title=element_blank())+
        #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
       guides(fill=guide_legend(override.aes=list(shape=c(21,21), 
            size=c(5,5), fill=c('brown3','blue3'))))
       #ggtitle("PC1 + PC3")
        #        shape=FALSE,
        #        size=FALSE)+
        #scale_size_manual(values=c(7,5))
#a2


pdf("~/oyster/figures/pca_scale.pdf", height=2.9, width=5)

ggarrange(a, a2, common.legend=TRUE, ncol=2, nrow=1, labels="AUTO")

dev.off()



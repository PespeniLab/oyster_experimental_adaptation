
require(corrplot) ; require(ape)
 #source the baypass R functions (check PATH)
 source("~/bin/baypass_2.1/utils/baypass_utils.R")

 #upload estimate of omega
 omega_33=as.matrix(read.table("~/oyster/analysis/baypass/33_core_mat_omega.out"))
 pop.names=c("LL", "OP", "TB")

 dimnames(omega_33)=list(pop.names,pop.names)

#Compute the correlation matrix
 cor.mat.33=cov2cor(omega_33)
tree.33=as.phylo(hclust(as.dist(1-cor.mat.33)))

cor.mat.all=cov2cor(omega_all)
tree.33=as.phylo(hclust(as.dist(1-cor.mat.33)))
tree2 <- rotate(tree.33, node=4)
tree3 <- rotate(tree2, node=5)

pdf("~/oyster/figures/hier_clust.pdf",height=3.35, width=3.35)
par(mfrow = c(1, 1), mar=c(2, 2, 2, 1), mgp=c(3, 1, 0), las=0)

plot(tree3,type="p", edge.width=2, cex=1.5)

dev.off()


 omega_all=as.matrix(read.table("~/oyster/analysis/baypass/all_core_mat_omega.out"))
 pop.names=c("LL_33","LL_5", "OP_33", "OP_5", "TB_33", "TB_5")
dimnames(omega_all)=list(pop.names,pop.names)

#Compute the correlation matrix
cor.mat.all=cov2cor(omega_all)

tree.all=as.phylo(hclust(as.dist(1-cor.mat.all)))

pdf("~/oyster/figures/hier_clust_allPops.pdf",height=3.35, width=3.35)
par(mfrow = c(1, 1), mar=c(2, 2, 2, 1), mgp=c(3, 1, 0), las=0)

plot(tree.all,type="p", edge.width=2, cex=1.5)
dev.off()



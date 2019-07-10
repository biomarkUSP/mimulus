library(qtl)
mimulus_qtl <- read.cross(format = "mm",file = "m_feb06.raw", mapfile="test_sum.txt")
mimulus_jm= jittermap(mimulus_qtl)
mimulus.qtl.probs = calc.genoprob(mimulus_qtl, step=1, error.prob = 0.001)

set.seed(1431)
mimulus_qtl_perm1= scanone(mimulus.qtl.probs, pheno.col = 4, method = "em", n.perm = 250)
set.seed(1432)
mimulus_qtl_perm2= scanone(mimulus.qtl.probs, pheno.col = 4, method = "em", n.perm = 250)
set.seed(1433)
mimulus_qtl_perm3= scanone(mimulus.qtl.probs, pheno.col = 4, method = "em", n.perm = 250)
set.seed(1434)
mimulus_qtl_perm4= scanone(mimulus.qtl.probs, pheno.col = 4, method = "em", n.perm = 250)

mimulus_qtl_perm = c(mimulus_qtl_perm1,mimulus_qtl_perm2,
                     mimulus_qtl_perm3,mimulus_qtl_perm4)
rm(mimulus_qtl_perm1,mimulus_qtl_perm2,
   mimulus_qtl_perm3,mimulus_qtl_perm4);

threshold_im = summary(mimulus_qtl_perm, alpha = 0.1)

{plot(mimulus_qtl_perm) 
  abline(v=threshold_im, lty= 2, col= "red")}



qtl_im= makeqtl(mimulus.qtl.probs, chr = sum_qtl_em$chr, pos = sum_qtl_em$pos, what = "prob")
plot(qtl_im)

mimulus_qtl.em= scanone(pheno.col = 4, mimulus.qtl.probs, method = "em")
sum_qtl_em = summary(mimulus_qtl.em, perm= mimulus_qtl_perm, alpha = 0.1, pvalues = TRUE)
sum_qtl_em


ggplot(data = mimulus_qtl.em,aes(x = pos, y = lod,colour=chr))+geom_line(size=1.3)+
  facet_wrap(~chr,nrow = 4,scales="free_x")+
  theme_minimal()+ ggtitle(label = "IM QTL model",
                           subtitle = "Leaf lenght (phenotype 4)")+
  ylab("LOD\n")+xlab("\nPoisition(cM)")+
  geom_hline(yintercept = threshold_im[1],linetype="dashed",size=.7)+
  theme(panel.border = element_rect(size=1.1, fill = NA,colour="black"),
        strip.text = element_text(face="bold"));

# analisando o cromossomo 12

mimulus_QTL= sim.geno(mimulus.qtl.probs, n.draws = 16, step = 1, off.end = 0, 
                     error.prob = 0.0001, map.function = "kosambi")

chr_12_effect= effectscan(mimulus_QTL,main="Chromossome 12",  pheno.col=4, chr = 12, get.se=FALSE, draw =TRUE, gap = 25,                                          mtick=c("line","triangle"),add.legend=TRUE, alternate.chrid=FALSE, 
                          ylab="QTL Effect", xlab="Linkage Group") 

chr_9_effect= effectscan(mimulus_QTL,main="Chromossome 9", pheno.col=4, chr = 9, get.se=FALSE, draw =TRUE, gap = 25,                                          mtick=c("line","triangle"),add.legend=TRUE, alternate.chrid=FALSE, 
                          ylab="QTL Effect", xlab="Linkage Group") 

summary(chr_12_effect)


### CIM

stepwiseqtl(mimulus.qtl.probs, pheno.col = 4, method = "hk", covar = NULL, model = "normal",
            incl.markers = TRUE, refine.locations = TRUE, additive.only = FALSE, scan.pairs = FALSE,
            keeplodprofile=TRUE, keeptrace=FALSE, verbose=TRUE, require.fullrank=FALSE)


name chr pos n.gen
Q1  9@49.0   9  49     3
Q2 12@13.0  12  13     3



#mimulus_cim_perm= cim(mimulus_cim_prob, method = "hk", n.marcovar = 3, pheno.col = 4, n.perm = 1000)

mimulus4.scan2 = scantwo(mimulus.qtl.probs,pheno.col = 4,method="hk",verbose=F)
save(mimulus4.scan2,file="mimulus4scan2.RData")
plot(mimulus4.scan2)
set.seed(1431)
mimulus_qtl_cim1= cim(mimulus.qtl.probs, pheno.col = 4, method = "hk", n.perm = 250)
set.seed(1432)
mimulus_qtl_cim2= cim(mimulus.qtl.probs, pheno.col = 4, method = "hk", n.perm = 250)
set.seed(1433)
mimulus_qtl_cim3= cim(mimulus.qtl.probs, pheno.col = 4, method = "hk", n.perm = 250)
set.seed(1434)
mimulus_qtl_cim4= cim(mimulus.qtl.probs, pheno.col = 4, method = "hk", n.perm = 250)

mimulus_qtl_cim = c(mimulus_qtl_cim1,mimulus_qtl_cim2,
                     mimulus_qtl_cim3,mimulus_qtl_cim4)
rm(mimulus_qtl_cim1,mimulus_qtl_cim2,
   mimulus_qtl_cim3,mimulus_qtl_cim4);

threshold_hk_cim =summary(mimulus_qtl_cim, alpha= c(0.1))
# 4.02

{plot(mimulus_qtl_cim) 
  abline(v=threshold_hk_cim, lty = 2, col = "red")}


mimulus_qtlcim= cim(mimulus.qtl.probs, method = "hk", pheno.col = 4, n.marcovar = 3, map.function = "kosambi")


sum_qtl_cim= summary(mimulus_qtl_cim, perms = mimulus_cim_perm, alpha = 0.1, pvalues = TRUE)
{plot(mimulus_qtlcim,col= "grey")
  abline(h = threshold_hk_cim,lty=2, col= "red")}

require(ggplot2)
ggplot(data = mimulus_qtlcim,aes(x = pos, y = lod,colour=chr))+geom_line(size=1.3)+
  facet_wrap(~chr,nrow = 4,scales="free_x")+
  theme_minimal()+ ggtitle(label = "CIM QTL model",
                           subtitle = "Leaf lenght (phenotype 4)")+
  ylab("LOD\n")+xlab("\nPoisition(cM)")+
   geom_hline(yintercept = threshold_hk_cim[1],linetype="dashed",size=.7)+
  theme(panel.border = element_rect(size=1.1, fill = NA,colour="black"),
        strip.text = element_text(face="bold"));


#plot(sum_qtl_cim)
qtl_cim = makeqtl(mimulus.qtl.probs, chr=c(9), pos=c(49), what="prob")


plot(qtl_cim, main="Composite Interval Mapping for Leaf Length (phenotype 4)");
plot(qtl_im, main="Interval Mapping for Leaf Length (phenotype 4)"); # 731 x 575

### MIM


out_fit1 <- fitqtl(mimulus.qtl.probs, pheno.col=16, qtl=qtl_im, method="hk", get.ests=TRUE)
summary(out_fit1)


effectplot(mimulus_QTL,mname1 = "9@49.0",mname2 = "12@13.0",pheno.col = 4)


## teste de modelos:

out_fit1 = fitqtl(mimulus.qtl.probs, pheno.col=16, qtl=qtl_im, method="hk", get.ests=TRUE,
                   formula = y~Q1+Q2);
sumamry(out_fit1)
out_fit2 = fitqtl(mimulus.qtl.probs, pheno.col=16, qtl=qtl_im, method="hk", get.ests=TRUE,
                   formula = y~Q1+Q2+Q1*Q2)

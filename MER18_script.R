#######################################################
# Cordier et al., 2018, Mol Ecol Res
#######################################################

# environment variables

# set path to the folder containing script and data
setwd("/path/to/the/folder")

# library dependencies
library(ranger)
library(vegan)
library(irr)
library(phyloseq)
library(metagenomeSeq)
## for parralel computation --> set the amount of cores available
library(doMC)
registerDoMC(cores = 8)

# for biotic indices calculation
library(devtools)
install_github("trtcrd/BBI")
library(BBI)

# function dependencies
source("sml_compo.R")
source("plot_ml.R")
source("bar_plot.R")

##### importing data 
# here we have to run the script for each OTU table
# note that the bacterial dataset has already been filtered, so does not need to be abundance filterd 
table <- "data/OTUtable_euksV1V2.txt"

comp <- read.table("data/metadata.txt", header=TRUE, sep="\t", dec=",")
otu <- read.table(table, header=TRUE, sep="\t", row.names=1)

s <- "samples_names"

## check if empty rows or columns at the end of file
otu[(dim(otu)[1]-10):(dim(otu)[1]),1:3]  ## normally there is nothing 
otu[1:3,(dim(otu)[2]-10):(dim(otu)[2])]  ## normally Taxa and X

### extract the OTU taxonomy assignments and total abundance from the OTU table as generated from slim
taxo <- otu[,c(1,dim(otu)[[2]]-1)]
dimnames(taxo)[[1]] <- dimnames(otu)[[1]]

### and keep the otu table only
otu <- data.frame(t(otu[,c(2:(dim(otu)[2]-2))]))
taxo<-  data.frame(t(taxo))

### subset to match the samples of OTU table (change the sample names accordingly)
comp <- subset(comp, comp[,s]!="")

## check the colnames and comp$samples_names setdiff?
#rownames(otu) == gsub("-",".",comp[,s])
comp <- subset(comp, gsub("-",".",comp[,s]) %in% rownames(otu))

### keep only samples with associated BIs values (if we have AMBI, we got the others)
otu  <- as.data.frame(subset(otu, is.na(comp$AMBI) == F))
comp <- subset(comp, is.na(comp$AMBI) == F)

# if things are right, we should have only TRUE
gsub("-",".",comp[,s]) == rownames(otu)

### keep the dataset at this point for plots and BIs calculation from metazoan assigned if eukaryotes marker
otu_raw <- otu
# plot sequencing depth
rsums<- rowSums(otu_raw)
plot(rsums)
## get samples with sequencing depth above the total average of reads per samples
seq_depth_cutoff <- 10000
otu  <- subset(otu, rsums >= seq_depth_cutoff)
comp <- subset(comp, rsums >= seq_depth_cutoff)

### Removing rares OTUs -- 100 for all markers but 1000 for V9 -- bacteria were already filtered!
otu <- otu[,colSums(otu)>100]
taxo <- taxo[,colnames(otu)]

### CSS normalization of the OTU table

### normalization with metagenomeSeq package
otu_to_norm  <- otu
comp_to_norm <- comp
samples_names <- comp[,s]

comp_ <- comp_to_norm[,c("Locality", "Station", "Grab")]
comp_$Station <- as.factor(comp_$Station)
comp_$Grab <- as.factor(comp_$Grab)
rownames(comp_) <- samples_names
rownames(otu_to_norm) <- samples_names

obj <- phyloseq(otu_table(otu_to_norm, taxa_are_rows = F), sample_data(comp_))
obj_m <- phyloseq_to_metagenomeSeq(obj)
p = cumNormStatFast(obj_m)
dat_mol_norm = cumNorm(obj_m, p = p)
otu_NORM <- t(MRcounts(dat_mol_norm, norm = TRUE, log = TRUE))

##### the data is ready!
OTU <- data.frame(otu_NORM)
OTU_raw <- otu
COMP <- comp

#### NMDS 

mds <-metaMDS(OTU, distance = "bray", k = 2, binary=F, trymax = 50, autotransform =F,    
              noshare = 0.1, wascores = T, expand = T, trace = 1,
              plot = F, old.wa = FALSE, display = "sites") 

quartz(width = 10, height=7)
par(mar=c(4, 4, 4, 12), xpd=TRUE)
plot(mds$points, col=as.numeric(COMP$col_plot), pch=COMP$col_plot, main = "NMDS on bray-curtis matrix")

env_dat <- COMP[,c("Distance_cage_gps", "AMBI", "NSI", "ISI", "NQI1")]
fit <- envfit(mds, env_dat, perm = 999, display = "sites", na.rm =T)
plot(fit, p.max = 0.05, col = "red")
ordisurf(mds, COMP[,"Distance_cage_gps"], add = TRUE, col="black")
leg <- unique(as.vector(COMP$Locality))
legend("topright", inset=c(-0.3,0), leg, col=unique(as.numeric(COMP$col_plot)), pch=unique(COMP$col_plot), box.lty=0)


##### SML

# which marker is it?
plot_title <- "V1V2 eukaryotes"

# how to agregate the data for taxonomic plot?
agregg <- c("Station", "Locality")

# this needs to be adjusted to the rank level to split the data
tax_rank <- 4

## get the taxonomic list at different taxo rank
tax_list <- bar_plot(OTU_raw, COMP, agregg, t(taxo), title_ = plot_title, font_size_ = 0.8, tax_rank = tax_rank, return_levels = T, remove_rare = T)
if (tax_list[length(tax_list)] == "zz_Unassigned")  tax_list[length(tax_list)] <- "Unassigned"
tax_list <- c(tax_list, "All")

## prepare table for gathering R2 and kappa by taxo group
statR <- array(NA, c(length(tax_list), 6))
statK <- array(NA, c(length(tax_list), 6))
dimnames(statR)[[2]] <- c("OTUs","Reads","AMBI", "ISI", "NSI", "NQI1")
dimnames(statR)[[1]] <- tax_list
dimnames(statK)[[2]] <- c("OTUs","Reads","AMBI", "ISI", "NSI", "NQI1")
dimnames(statK)[[1]] <- tax_list

for (k in tax_list)
{
  if (k == "All") 
  {
    tt <- t(taxo)
    OTUtpML <- OTU_raw
    OTUtpPL <- OTU
    TAXOtp <- tt
    table_ <- cbind(TAXOtp, t(OTUtpPL))
    
  } else {
    tt <- t(taxo)
    ### grep need to be at the rank we split... 
    #tp <- c()
    #for (i in 1:length(tt[,2])) tp <- c(tp, unlist(strsplit(as.character(tt[i,"Taxa"]), split=";", fixed=TRUE))[tax_rank_for_list])
    tmp <- grep(k, tt[,2], ignore.case=F)
    #tmp <- grep(k, tmp)
    OTUtpML <- OTU[,tmp]
    OTUtpPL <- OTU_raw[,tmp]
    TAXOtp <- tt[tmp,]
    ### export the table # length(subset(colnames(OTUtpPL) == rownames(TAXOtp), (colnames(OTUtpPL) == rownames(TAXOtp)) == T))
    table_ <- cbind(TAXOtp, t(OTUtpPL))
    
  }
  # for dataset that will contain only one OTUs...
  if (is.null(dim(OTUtpML)) == F)
  {
    bar_plot(OTUtpPL, COMP, agregg, TAXOtp, title_ = paste(k,plot_title, "-", ncol(OTUtpML),
                                                           "OTUs","-",sum(OTUtpPL), "reads"), font_size_ = 0.8, tax_rank = tax_rank, lastrank=F, pdf =T, taxo_group=k, remove_rare = F)
    bar_plot(OTUtpPL, COMP, agregg, TAXOtp, title_ = paste(k,plot_title, "-", ncol(OTUtpML),
                                                           "OTUs","-",sum(OTUtpPL), "reads"), font_size_ = 0.8, tax_rank = tax_rank, lastrank=T, pdf =T, taxo_group=k, remove_rare = F)
    statR[k,"OTUs"] <- statK[k,"OTUs"] <- ncol(OTUtpML)
    statR[k,"Reads"] <- statK[k,"Reads"] <- sum(OTUtpPL)
  } else {
    statR[k,"OTUs"] <- statK[k,"OTUs"] <- 1
    statR[k,"Reads"] <- statK[k,"Reads"] <- sum(OTUtpPL)
  }
  for (index in c("AMBI", "ISI", "NSI", "NQI1"))
  {
    test <- sml_compo(OTUtpML, COMP, index)
    test_c <- plot_ml(test, COMP, index, paste("RF - ", plot_title, k), pdf=T, k)
    statR[k,index] <- test_c$R2
    statK[k,index] <- test_c$KAP
    
    ## plot importance
    if (is.null(dim(OTUtpML)) == F) # for dataset that will contain only one OTUs...
    {
      mod <- ranger(COMP[,index] ~ ., data=OTUtpML, mtry=floor(dim(OTUtpML)[2]/3), num.trees = 300, importance= "impurity", write.forest = T)
      imp <- tail(sort(mod$variable.importance), 10)
      tx <- TAXOtp[c(names(imp)), "Taxa"]
      for (i in 1:length(tx)) tx[i] <- tail(unlist(strsplit(as.character(tx[i]), split=";", fixed=TRUE)), 1)
      #quartz(width = 5, height=4)
      pdf(width = 5, height=4, file = paste(k, "/", paste(index, "_",k,"_import_plot.pdf", sep=""), sep=""))
      mat <- matrix(c(2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,11,1,12,1,13,1,14,1,15,1), nrow=14, byrow = T)
      layout(mat, widths=c(1,9,1,1,1,1,1,1,1,1,1,1,1,1,1,1), heights=c(2.4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
      #layout.show(n=15)
      p <- barplot(imp, horiz = T, xlab="Variable importance", main=paste(k, "OTUs importance for", index), las=2, cex.names = 0.7)
      text(0, p, labels=tx, pos=4, cex=0.7)
      mtext(colSums(OTUtpPL[,c(names(tx))]), side=2, at=p-0.4, cex=0.3, las=1, line=1)
      par(mar=c(0.2,0.8,0.2,0.8), xpd=TRUE)
      plot(0,type="n", axes=F, xlab="", ylab="")
      for (i in names(rev(imp))) plot(OTU[,i], COMP[,index], cex=0.3, pch=16, yaxt="n", xaxt="n")
      dev.off()
    }
  }
  ## export the OTU table for each sub-groups
  write.table(table_, paste(k, "/", k, "_otu_table_with_assign.txt", sep=""), col.names = T, row.names = T, quote=F, sep="\t")
}

# export the summarized stats
write.table(statR, paste(plot_title, "_R2.txt", sep=""), col.names = T, row.names = T, quote=F, sep="\t")
write.table(statK, paste(plot_title, "_KAP.txt", sep=""), col.names = T, row.names = T, quote=F, sep="\t")


### Compute BI indices with BBI package

# first fetch the metazoan OTUs
tmp <- grep("Metazoa", t(taxo)[,2])
OTU_metaz <- OTU_raw[,tmp]
tax <- t(taxo[,tmp])

# concatenate to make the data.frame for BBI
metaz_data <- cbind(Taxa = tax[,"Taxa"], t(OTU_metaz))

# calculation of BIs values
BIs <- BBI(metaz_data)

# plot the results with the reference BIs values obtianed from macrofauna
eDNA_AMBI <- plot_ml(BIs$BBI[,"AMBI"], COMP, "AMBI", paste("Taxonomy-based - ", plot_title))
eDNA_NSI <- plot_ml(BIs$BBI[,"NSI"], COMP, "NSI", paste("Taxonomy-based - ", plot_title))
eDNA_ISI <- plot_ml(BIs$BBI[,"ISI"], COMP, "ISI", paste("Taxonomy-based - ", plot_title))
eDNA_NQI1 <- plot_ml(BIs$BBI[,"NQI1"], COMP, "NQI1", paste("Taxonomy-based - ", plot_title))



##########################################################





#### barplot of taxonomic rank relative abundance

bar_plot <- function(otu_table, comp_, aggreg = F, taxo_, title_, font_size_ = 0.8, horizontal_ = T, tax_rank, metazoan = F, lastrank = F , pdf = F, taxo_group = NULL, return_levels = F, remove_rare = T, remove_unass = F) {

  otu <- otu_table
  comp <- comp_
  tax <- taxo_
  title <- title_
  font_size <- font_size_
  horizontal <- horizontal_
  tax_rank <- tax_rank
  metazoan <- metazoan
  lastrank <- lastrank
  pdf <- pdf
  taxo_group <- taxo_group
  
  if (lastrank) { taxRank <- "last" } else { taxRank <- tax_rank }
  
  ## title
  if (lastrank) title <- paste(title, "(last rank)")
  if (!lastrank) title <- paste(title, " (rank ", tax_rank, ")", sep="")
  # 
  # otu <- OTUtpPL
  # comp <- COMP
  # tax <- TAXOtp
  # title <- paste(k,taxo_group, "-", ncol(OTUtpML),"OTUs","-",sum(OTUtpPL), "reads")
  # font_size <- 0.8
  # horizontal_ <- T
  # tax_rank <- tax_rank
  # metazoan <- F
  # lastrank <- F
  # pdf <- T
  # taxo_group <- k

  # need to remove [ ] or weird character in taxo file if any
  for (i in 1:length(tax)) tax[i] <- gsub("[^[:alnum:][:blank:]_;+?&/\\-]", "", tax[i], c) 

  # tax_rank 2 if forams..
  # tax rank 7 if metazoa 
  # tax rank 2 if bacteria
  # tax rank 4 if eukaryotes
  
  ## if metazoan keep only metazoan
  if (metazoan) 
  {
    metaz <- grep("Metazoa", tax[,"Taxa"])
    otu <- otu[,metaz]
    tax <-tax[metaz,]
    tax_rank <- 7
  }
  
  ## if remove unassigned
  if (remove_unass) 
  {
    unass <- grep("Unassigned", tax[,"Taxa"])
    otu <- otu[,-unass]
    tax <-tax[-unass,]
    #tax_rank <- 7
  }
  
  # merging the data by grab 
  #otu_ <- aggregate(otu, by=list(comp, comp$Station, comp$Locality), FUN=sum)
  #otu_ <- aggregate(otu, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=sum)
  
  if (aggreg==F) dat <- otu
  if (aggreg!=F)
  {
    if (length(aggreg)==1) otu_ <- aggregate(otu, by=list(comp[,aggreg[1]]), FUN=sum)
    if (length(aggreg)==2) otu_ <- aggregate(otu, by=list(comp[,aggreg[1]], comp[,aggreg[2]]), FUN=sum)
    if (length(aggreg)==3) otu_ <- aggregate(otu, by=list(comp[,aggreg[1]], comp[,aggreg[2]], comp[,aggreg[3]]), FUN=sum)
    # and reformatting to get single field with locality > station > grab
    dat <- otu_[,4:dim(otu_)[[2]]]
  }
  ### then sorting by station and locality 
  #dat <- otu_[with(otu_, order(Group.2, Group.3)),]
  
  # and repasting the names with all infos
  ## if reordering by station and locality ## rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[with(otu_, order(Group.2, Group.3)),c("Group.3","Group.2","Group.1")]), stringsAsFactors=FALSE))
  if (aggreg!=F)
  {
    if (length(aggreg)==1) rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[,c("Group.1")]), stringsAsFactors=FALSE))
    if (length(aggreg)==2) rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[,c("Group.2","Group.1")]), stringsAsFactors=FALSE))
    if (length(aggreg)==3) rownames(dat) <- do.call(paste, as.data.frame(as.matrix(otu_[,c("Group.3","Group.2","Group.1")]), stringsAsFactors=FALSE))
  }
  # extract the "family" assignment
  if (lastrank)
  {
    taxa <- seq(from = 1, to = dim(dat)[2])
    for (i in 1:dim(dat)[2])
    {
      if (tax[i,"Taxa"] == "") tax[i,"Taxa"] <- "Unassigned"
      tmp <- tail(unlist(strsplit(as.character(tax[i,"Taxa"]), split=";", fixed=TRUE)), 1)
      taxa[i] <- tmp
      if (tmp == "Unassigned") taxa[i] <- "zz_Unassigned"
    }
  } else
  {
    taxa <- seq(from = 1, to = dim(dat)[2])
    for (i in 1:dim(dat)[2])
    {
      ### if the taxonomic assignment is empty -> Unassigned
      if (tax[i,"Taxa"] == "") tax[i,"Taxa"] <- "Unassigned" 
      ### extract the tax_rank of the assignment
      tmp <- unlist(strsplit(as.character(tax[i,"Taxa"]), split=";", fixed=TRUE))[tax_rank] 
      ### if it is not empty
      if (is.na(tmp) == FALSE)
      {
        taxa[i] <- tmp
        if (tmp == "Unassigned") taxa[i] <- "zz_Unassigned"
      } else {
        ### extract only the tax_rank, if lower, it goes to unassigned...
        # tmp <- tail(unlist(strsplit(as.character(tax[i,"Taxa"]), split=";", fixed=TRUE)), 1) ### previous code
        #tmp <- unlist(strsplit(as.character(tax[i,"Taxa"]), split=";", fixed=TRUE))[tax_rank]
        #if (is.na(tmp) == TRUE) taxa[i] <- "zz_Unassigned"
        #if (tmp == "Unassigned") taxa[i] <- "zz_Unassigned"
        taxa[i] <- "zz_Unassigned"
      }
    }
  }
  # creating a table for sum family
  dat_ <- as.data.frame(array(NA, c(dim(dat)[1], length(levels(as.factor(taxa))))))
  dimnames(dat_)[[1]] <- dimnames(dat)[[1]]
  ll <- levels(as.factor(taxa))
  for (i in 1:length(ll)) ll[i] <- strsplit(ll[i], split="(", fixed=TRUE)[[1]][1]    # need to remove parenthesis if any
  dimnames(dat_)[[2]] <- ll
  
  # summing the abundance by tax_rank
  for (i in ll)
  {
    tt  <- dat[,grep(i, taxa)]
    if (is.data.frame(tt)) { dat_[,i] <- rowSums(tt) } else { dat_[,i] <- tt }
  }
  
  # normalizing the otu table for % 
  dat_ <- dat_ / rowSums(dat_) *100
  dat_[dat_ == "NaN"] <- 0
  
  ## remove rare families
  if (remove_rare) 
  {
    if (ncol(dat_) > 1) dat_ <-  dat_[,colSums(dat_)>1]
    # REnormalizing the otu table for % 
    dat_ <- dat_ / rowSums(dat_) *100
    dat_[dat_ == "NaN"] <- 0
  }
  
  
  #barplot
  col <- colorRampPalette(c("black", "yellow", "blue", "darkorchid1", "red", "green", "orange", "white", "lightblue1", "lavenderblush", "olivedrab4", "seagreen1", "grey"), bias=1, interpolate = "linear")(dim(dat_)[2])
  
  #col <- colorRampPalette("palette")(dim(dat_)[2])
  
  if (horizontal_) 
  {
    if (pdf)
    {
      if (is.null(taxo_group)) 
      {
        pdf(file = "taxo_group_plot.pdf")
      } else {
        dir.create(file.path(getwd(), taxo_group))
        pdf(file = paste(taxo_group, "/", paste(taxo_group, "_rank_", taxRank, "_plot.pdf", sep=""), sep=""))
      }
      par(mar=c(4, 8, 2, 10), xpd=TRUE)
      barplot(as.matrix(t(dat_)), col=col, main=title, las=1, horiz=TRUE, cex.names = font_size, border = TRUE, xlab="relative abundance (%)")
      legg <- dimnames(dat_)[[2]]
      if(!remove_unass) legg[length(legg)] <- "Unassigned"
      legend("topright",legend = legg,fill = col, cex=0.6, inset=c(-0.55,0.02), box.lty=0)
      dev.off()
    } else
    {
      quartz()
      par(mar=c(4, 8, 2, 10), xpd=TRUE)
      barplot(as.matrix(t(dat_)), col=col, main=title, las=1, horiz=TRUE, cex.names = font_size, border = TRUE, xlab="relative abundance (%)")
      legg <- dimnames(dat_)[[2]]
      if(!remove_unass) legg[length(legg)] <- "Unassigned"
      legend("topright",legend = legg,fill = col, cex=0.6, inset=c(-0.55,0.02), box.lty=0)
    }
  } else {
    if (pdf)
    {
      if (is.null(taxo_group)) 
      {
        pdf(file = "taxo_group_plot.pdf")
      } else {
        dir.create(file.path(getwd(), taxo_group))
        pdf(file = paste(taxo_group, "/", paste(taxo_group, "_rank_", taxRank, "_plot.pdf", sep=""), sep=""))
      }
      par(mar=c(8, 4, 2, 10), xpd=TRUE)
      barplot(as.matrix(t(dat_)), col=col, main=title, las=3, horiz=FALSE, cex.names = font_size, border = TRUE, ylab="relative abundance (%)")
      legg <- dimnames(dat_)[[2]]
      if(!remove_unass) legg[length(legg)] <- "Unassigned"
      legend("topright",legend = legg,fill = col, cex=0.6, inset=c(-0.35,0), box.lty=0)
      dev.off()
    } else
    {
      quartz()
      par(mar=c(8, 4, 2, 10), xpd=TRUE)
      barplot(as.matrix(t(dat_)), col=col, main=title, las=3, horiz=FALSE, cex.names = font_size, border = TRUE, ylab="relative abundance (%)")
      legg <- dimnames(dat_)[[2]]
      if(!remove_unass) legg[length(legg)] <- "Unassigned"
      legend("topright",legend = legg,fill = col, cex=0.6, inset=c(-0.35,0), box.lty=0)
    }
  }
  
  if (return_levels) return(levels(as.factor(dimnames(dat_)[[2]])))
  
}

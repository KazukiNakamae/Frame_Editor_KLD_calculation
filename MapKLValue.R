#-- Load library ----------------------------------------------------#
message("---Load libraries---")

library(CrispRVariants)
library(Biostrings)
library(ggalluvial)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(exactRankTests)


classifier.res.dir <- "./summary"
aligner.res.dir <- "./variant_call_result"
output.prefix <- "./MapKLValue_output"
condition.label.vec <- c("WT1", "WT2", "WT3", "T41", "T42", "T43", "GB1", "GB2", "GB3", "HifiH1", "HifiH2", "HifiH3", "HifiH41", "HifiH42", "HifiH43", "HifiHB1", "HifiHB2", "HifiHB3")
ignore.list <- ""
indelphi.file.path <- "/Volumes/databank4/nakade_mutatiobn_analysis/20211102/NGS_WB_W4/Indelphi_files"


#-- Make output directory ----------------------------------------------------#
message("---Make directory for the analysis---")

time.id <- format(Sys.time(), "%Y%m%d%H%M%S")
output.name <- paste("MapKLValue_at_", time.id, sep = "")
dir.create(file.path(output.prefix), showWarnings=FALSE)
output.dir <- file.path(output.prefix, output.name)
dir.create(file.path(output.dir), showWarnings=FALSE)

################################################################################################################################
### result file information data frame
################################################################################################################################
classifier.prefix.path.vec <- list.dirs(path = file.path(classifier.res.dir), full.names = TRUE, recursive = FALSE)
classifier.res.path.vec <- character(0)
for(classifier.prefix.path in classifier.prefix.path.vec){
    classifier.res.path.vec <- c(classifier.res.path.vec, list.dirs(path = file.path(classifier.prefix.path), full.names = TRUE, recursive = FALSE))
}
classifier.res.sample.name.vec <- sapply(classifier.res.path.vec, function(x){
    paste(c(tail(strsplit(x, "_")[[1]], n=1)[1], paste(tail(strsplit(tail(strsplit(x, "/")[[1]], n=2)[1], "_")[[1]], n=2), collapse="")), collapse="-")
    })
classifier.res.sample.target.vec <- sapply(classifier.res.sample.name.vec, function(x){head(strsplit(x, "-")[[1]], n=1)})
classifier.res.sample.label.vec <- sapply(classifier.res.sample.name.vec, function(x){tail(strsplit(x, "-")[[1]], n=1)})
classifier.res.file.df <- data.frame(
  classifier.res.dir.path = classifier.res.path.vec
  , sample.name = classifier.res.sample.name.vec
  , sample.target = classifier.res.sample.target.vec
  , sample.label = classifier.res.sample.label.vec
  , stringsAsFactors = FALSE
)

aligner.prefix.path.vec <- list.dirs(path = file.path(aligner.res.dir), full.names = TRUE, recursive = FALSE)
aligner.res.path.vec <- character(0)
for(aligner.prefix.path in aligner.prefix.path.vec){
    aligner.res.path.vec <- c(aligner.res.path.vec, list.dirs(path = file.path(aligner.prefix.path), full.names = TRUE, recursive = FALSE))
}
aligner.res.sample.name.vec <- sapply(aligner.res.path.vec, function(x){
    paste(c(tail(strsplit(x, "_")[[1]], n=1)[1], paste(tail(strsplit(tail(strsplit(x, "/")[[1]], n=2)[1], "_")[[1]], n=2), collapse="")), collapse="-")
    })
aligner.res.file.df <- data.frame(
  aligner.res.dir.path = aligner.res.path.vec
  , sample.name = aligner.res.sample.name.vec
  , stringsAsFactors = FALSE
)
res.file.df <- merge(classifier.res.file.df, aligner.res.file.df, by = "sample.name", all = TRUE)
################################################################################################################################
### Filtering
################################################################################################################################
message("Filtering data")
if(ignore.list != ""){
  remove.target.vec <- unlist(read.csv(ignore.list, header = FALSE, stringsAsFactors = FALSE)[1,])
}else{
  remove.target.vec <- character(0)
}

# Remove gene contains low number of reads
message("Check amount of reads...")
has.sufficint.reads.vec <- apply(res.file.df, MARGIN = 1, function(row){
  all.df <- read.csv(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"), header = TRUE)

  if(sum(all.df$X.Reads) < 1000){
    message(paste0(row["sample.name"], " dosen't have sufficint reads (<1,000 reads)."))
    return(FALSE)
  }else{
    return(TRUE)
  }
})
remove.target.vec <- c(remove.target.vec, unique(res.file.df$sample.target[!has.sufficint.reads.vec]))

remove.target.vec <- unique(remove.target.vec)
message("Remove...")
print(remove.target.vec)
fltr.res.file.df <- filter(res.file.df, !(sample.target %in% remove.target.vec))
message("Filtering is done.")

# Read functions for calculating KL values and mapping them
script.dir <- "."
source(file.path(script.dir, "MapKLValue_functions.R"))

################################################################################################################################
message("---Run optional analysis---")
########Comparison between InDelphi InDel Variants and detected InDel Variants
message("Comparison between InDelphi InDel Variants and detected InDel Variants")
indelphi.data.list <- apply(res.file.df, MARGIN = 1, function(row){
  temp.file.name <- paste0("inDelphi_output_", as.character(row["sample.target"]), ".csv")
  if((row["sample.label"] %in% c(condition.label.vec[1])) & file.exists(file.path(indelphi.file.path, temp.file.name))){
    temp.data <- read.csv(file.path(indelphi.file.path, temp.file.name))
    temp.data$Category <- gsub("del", "D", temp.data$Category)
    temp.data$Category <- gsub("ins", "I", temp.data$Category)
    # add position information
    temp.data["position"] <- temp.data$Genotype.position - temp.data$Length - temp.data$Microhomology.length
    # if variant that does not have mh deletion, the value is 0 
    temp.data$Microhomology.length <- apply(temp.data, MARGIN = 1, function(row){
      if(is.na(row["Microhomology.length"])){
        row["Microhomology.length"] <- 0 # the insertion position of ambigous because we cannot indentify the insertion base from the duplication sequences
      }
      return(as.numeric(row["Microhomology.length"]))
    })
    # define position
    temp.data$position <- apply(temp.data, MARGIN = 1, function(row){
      if(is.na(row["position"])){
        row["position"] <- "" # the insertion position of ambigous because we cannot indentify the insertion base from the duplication sequences
      }else if(as.numeric(row["position"]) > -1){
        row["position"] <- as.numeric(row["position"]) + 1 # if the "0" position is nothing, add +1 in the positive value (including 0) case.
      }else{
        if(as.numeric(row["Microhomology.length"]) == 0){ # add -1 in case of the value is negative and nhej-mediated deletion.
          row["position"] <- as.numeric(row["position"]) - 1
        }
      }
      return(row["position"])
    })
    # annotate mmej mediated deletion and indertion
    mut.lebel <- paste0(temp.data$position, ":", temp.data$Length, temp.data$Category)
    for(mut.ind in 1:length(mut.lebel)){
      if(temp.data$Microhomology.length[mut.ind]>0){
        mut.lebel[mut.ind] <- paste0(mut.lebel[mut.ind], ":mh", as.character(temp.data$Microhomology.length[mut.ind]))
      }
      if(temp.data$Inserted.Bases[mut.ind]>0){
        mut.lebel[mut.ind] <- paste0(mut.lebel[mut.ind], ":", temp.data$Inserted.Bases[mut.ind])
      }
      mut.lebel[mut.ind] <- gsub("[[:blank:]]", "", mut.lebel[mut.ind])
    }
    # browser()

    # aggregated.freq.table <- aggregate(reads ~ label, data = aggregated.freq.table, sum)
    aggregated.freq.table <- data.frame(
      label = mut.lebel,
      genotype = temp.data$Genotype,
      reads = temp.data$Predicted.frequency
    )
    temp.df <- matrix(
      c(aggregated.freq.table$reads,
      aggregated.freq.table$genotype),
      ncol = 2
    )
    rownames(temp.df) <- aggregated.freq.table$label
    return(temp.df)
  }else if(!(row["sample.label"] %in% c(condition.label.vec[1]))){
    return(NULL)
  }else{
    message(paste0("Ignore ", row["sample.name"], " because ", temp.file.name, " was not found."))
    return(NULL)
  }
})
indelphi.data.name.vec <- res.file.df$sample.name[!sapply(indelphi.data.list, is.null)]
indelphi.data.list <- Filter(Negate(is.null), indelphi.data.list)
names(indelphi.data.list) <- indelphi.data.name.vec

genename.vec <- c("A2LD1","AR","ATP5G1","BRP44","C16orf93","C17orf96","C1orf204","CD34","FAM159B","HBXIP","MAP3K1","ODZ3","PCSK91","PPP1R12C-Exon2","PTPLB","RP4-697K147","TMEM57","TPH22","XRCC6","ZNF498")
total.diagonal.vec = numeric(0)
# Map HeatMap
for(condition.label in condition.label.vec){
  # browser()
  MakeHeatMap(
    fltr.res.file.df,
    indelphi.data.list,
    "Predicted InDel Variants based on InDelphi",
    paste0("Detected InDel Variants in ", condition.label, " sample (CRISPResso2 Classification)"),
    file.path(output.dir, paste0("Comparison_between_InDelphi_and_CRISPResso2_InDel_class_in_", condition.label, "_samples")),
    condition.label,
    class_type = "base"
  )
  # get diagonal element
  KL.mat <- as.matrix(read.table(paste0(file.path(output.dir, paste0("Comparison_between_InDelphi_and_CRISPResso2_InDel_class_in_", condition.label, "_samples")), ".csv"), sep=","))
  diagonal.vec <- KL.mat[col(KL.mat)==row(KL.mat)]
  names(diagonal.vec) <- rownames(KL.mat)
  for(genename in genename.vec){
    total.diagonal.vec <- c(total.diagonal.vec, diagonal.vec[genename])
  }
}

# make diagonal matrix 
total.diagonal.mat <- matrix(total.diagonal.vec, ncol=length(condition.label.vec))
colnames(total.diagonal.mat) <- condition.label.vec
rownames(total.diagonal.mat) <- genename.vec

# plot heatmap for diagonal KL value
for(gene.ind in 1:length(genename.vec)){
  min.diagonal <- min(total.diagonal.mat[gene.ind,][!is.na(total.diagonal.mat[gene.ind,])])
  max.diagonal <- max(total.diagonal.mat[gene.ind,][!is.na(total.diagonal.mat[gene.ind,])])
  data <- matrix(c(total.diagonal.mat[gene.ind,], total.diagonal.mat[gene.ind,]), nrow = length(total.diagonal.mat[gene.ind,]))
  rownames(data) <- condition.label.vec
  colnames(data) <- c("", "")
  plot.diagonal.heatmap(data, 2, 
    min.diagonal, max.diagonal, 
    file.path(output.dir, paste0("Diagonal_Comparison_between_InDelphi_and_CRISPResso2_InDel_class_in_", genename.vec[gene.ind])), 
    gene.name=genename.vec[gene.ind], 
    "Symmetrized KL divergence")
}

# WT
wt.median.diagonal.vec <- apply(total.diagonal.mat[, grepl("WT", colnames(total.diagonal.mat))], MARGIN = 1, median)
# T4
t4.median.diagonal.vec <- apply(total.diagonal.mat[, grepl("T4", colnames(total.diagonal.mat))], MARGIN = 1, median)
# GB
gb.median.diagonal.vec <- apply(total.diagonal.mat[, grepl("GB", colnames(total.diagonal.mat))], MARGIN = 1, median)
# HifiH
hifih.median.diagonal.vec <- apply(total.diagonal.mat[, grepl("HifiH[0-9]$", colnames(total.diagonal.mat))], MARGIN = 1, median)
# HifiH4
hifih4.median.diagonal.vec <- apply(total.diagonal.mat[, grepl("HifiH4", colnames(total.diagonal.mat))], MARGIN = 1, median)
# HifiHB
hifihb.median.diagonal.vec <- apply(total.diagonal.mat[, grepl("HifiHB", colnames(total.diagonal.mat))], MARGIN = 1, median)

new.condition.label.vec <- c("Cas9", "dFE", "iFE", "HiFi Cas9", "HiFi dFE", "HiFi iFE")

median.diagonal.mat <- matrix(c(wt.median.diagonal.vec, t4.median.diagonal.vec, gb.median.diagonal.vec, hifih.median.diagonal.vec, hifih4.median.diagonal.vec, hifihb.median.diagonal.vec), nrow=length(genename.vec))
colnames(median.diagonal.mat) <- new.condition.label.vec
rownames(median.diagonal.mat) <- genename.vec

# plot heatmap for diagonal KL value
for(gene.ind in 1:length(genename.vec)){
  min.diagonal <- min(median.diagonal.mat[gene.ind,][!is.na(median.diagonal.mat[gene.ind,])])
  max.diagonal <- max(median.diagonal.mat[gene.ind,][!is.na(median.diagonal.mat[gene.ind,])])
  data <- matrix(c(median.diagonal.mat[gene.ind,], median.diagonal.mat[gene.ind,]), nrow = length(median.diagonal.mat[gene.ind,]))
  rownames(data) <- new.condition.label.vec
  colnames(data) <- c("", "")
  # we can adjust scale so that the the value of Cas9 should be center of the scale.
  delta.min <- abs(min.diagonal - median.diagonal.mat[gene.ind,"Cas9"])
  delta.max <- abs(max.diagonal - median.diagonal.mat[gene.ind,"Cas9"])
  if(delta.min > delta.max){
    min.lim <- min.diagonal
    max.lim <- min.diagonal + ( 2 * delta.min )
  }else{
    min.lim <- max.diagonal - ( 2 * delta.max )
    max.lim <- max.diagonal
  }
  # plot
  plot.diagonal.heatmap(data, 2, 
    min.lim, 
    max.lim, 
    file.path(output.dir, paste0("Median_Diagonal_Comparison_between_InDelphi_and_CRISPResso2_InDel_class_in_", genename.vec[gene.ind])),
    gene.name=genename.vec[gene.ind], 
    "Median Symmetrized KL divergence")
}

# plot heatmap for diagonal KL value(Scale is common among genes)
min.diagonal <- min(median.diagonal.mat[!is.na(median.diagonal.mat)])
max.diagonal <- max(median.diagonal.mat[!is.na(median.diagonal.mat)])
for(gene.ind in 1:length(genename.vec)){
  data <- matrix(c(median.diagonal.mat[gene.ind,], median.diagonal.mat[gene.ind,]), nrow = length(median.diagonal.mat[gene.ind,]))
  rownames(data) <- new.condition.label.vec
  colnames(data) <- c("", "")
  # plot
  plot.diagonal.heatmap(data, 2, 
    min.diagonal, 
    max.diagonal, 
    file.path(output.dir, paste0("Scalefixed_Median_Diagonal_Comparison_between_InDelphi_and_CRISPResso2_InDel_class_in_", genename.vec[gene.ind])), 
    gene.name=genename.vec[gene.ind], 
    "Median Symmetrized KL divergence")
}

total.diagonal.mat.melt <- melt(total.diagonal.mat)

total.diagonal.mat.melt.p <- ggplot(total.diagonal.mat.melt, aes(x = Var2, y = value)) + labs(x = "Condition", y = "Symmetrized KL Divergence") + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.2)

ggsave(file = file.path(output.dir, "All_Diagonal_Comparison_between_InDelphi_and_CRISPResso2_InDel_class_boxplot.png"), plot = total.diagonal.mat.melt.p, dpi = 350, width = 16, height = 2)

t.test(total.diagonal.mat[, "WT1"], total.diagonal.mat[, "T41"], var.equal=F)
# p-value = 0.3726

t.test(c(total.diagonal.mat[, "WT1"], total.diagonal.mat[, "WT2"], total.diagonal.mat[, "WT3"]), c(total.diagonal.mat[, "T41"], total.diagonal.mat[, "T42"], total.diagonal.mat[, "T43"]), var.equal=F)
# p-value = 0.1695


# t-test of KL divergence in WT vs FE samples
res.t.test.df <- apply(total.diagonal.mat, MARGIN = 1, function(row){
  sample.regular.exp.vec <- c("T4", "GB", "HifiH[0-9]$", "HifiH4", "HifiHB")
  p.value.vec <- numeric(0)
  for(sample.regular.exp.ind in 1:length(sample.regular.exp.vec)){
    not.na.vec <- !is.na(row[grepl(sample.regular.exp.vec[sample.regular.exp.ind], names(row))])
    if(length(not.na.vec[not.na.vec == TRUE]) >2 ){
      res.t.test <- t.test(row[grepl("WT", names(row))], row[grepl(sample.regular.exp.vec[sample.regular.exp.ind], names(row))], var.equal=F)
      p.value.vec <- c(p.value.vec, res.t.test$p.value)
      if(res.t.test$p.value < 0.05){
        p.value.vec <- c(p.value.vec, TRUE)
      }else{
        p.value.vec <- c(p.value.vec, FALSE)
      }
    }else{
      p.value.vec <- c(p.value.vec, NA)
      p.value.vec <- c(p.value.vec, FALSE)
    }
  }
  return(p.value.vec)
})
row.names(res.t.test.df) <- c("dFE", "dFE_alpha_five", "iFE", "iFE_alpha_five", "HiFi Cas9", "HiFi Cas9_alpha_five", "HiFi dFE", "HiFi dFE_alpha_five", "HiFi iFE", "HiFi iFE_alpha_five")
res.t.test.df <- t(res.t.test.df)
SaveTable(res.t.test.df, file.path(output.dir, "All_Diagonal_Comparison_between_InDelphi_and_CRISPResso2_InDel_class_t_test_by_gene"))

# Exact Wilcoxon rank sum test of KL divergence in WT vs FE samples
res.wilcox.exact.df <- apply(total.diagonal.mat, MARGIN = 1, function(row){
  sample.regular.exp.vec <- c("T4", "GB", "HifiH[0-9]$", "HifiH4", "HifiHB")
  p.value.vec <- numeric(0)
  for(sample.regular.exp.ind in 1:length(sample.regular.exp.vec)){
    not.na.vec <- !is.na(row[grepl(sample.regular.exp.vec[sample.regular.exp.ind], names(row))])
    if(length(not.na.vec[not.na.vec == TRUE]) >2 ){
      res.wilcox.exact <- wilcox.exact(row[grepl("WT", names(row))], row[grepl(sample.regular.exp.vec[sample.regular.exp.ind], names(row))], paired=F)
      p.value.vec <- c(p.value.vec, res.wilcox.exact$p.value)
      if(res.wilcox.exact$p.value < 0.05){
        p.value.vec <- c(p.value.vec, TRUE)
      }else{
        p.value.vec <- c(p.value.vec, FALSE)
      }
    }else{
      p.value.vec <- c(p.value.vec, NA)
      p.value.vec <- c(p.value.vec, FALSE)
    }
  }
  return(p.value.vec)
})
row.names(res.wilcox.exact.df) <- c("dFE", "dFE_alpha_five", "iFE", "iFE_alpha_five", "HiFi Cas9", "HiFi Cas9_alpha_five", "HiFi dFE", "HiFi dFE_alpha_five", "HiFi iFE", "HiFi iFE_alpha_five")
res.wilcox.exact.df <- t(res.wilcox.exact.df)
SaveTable(res.wilcox.exact.df, file.path(output.dir, "All_Diagonal_Comparison_between_InDelphi_and_CRISPResso2_InDel_class_wilcox_exact_test_by_gene"))

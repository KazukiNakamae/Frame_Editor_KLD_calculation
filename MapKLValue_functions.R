SaveTable <- function(table, file.name.without.ext){
  options(warn = -1)
  if(is.list(table) == TRUE){
    if(is.null(unlist(table))){
      return()
    }
    write.table(data.frame(matrix(unlist(table), nrow=length(table), byrow=T)) , file = paste0(file.name.without.ext, ".csv"), quote=FALSE, col.names=TRUE, row.names=TRUE,append=TRUE, sep = ",")
  }else{
    write.table(table , file = paste0(file.name.without.ext, ".csv"), quote=FALSE, col.names=TRUE, row.names=TRUE,append=TRUE, sep = ",")
  }
  options(warn = 0)
  saveRDS(table , file = paste0(file.name.without.ext, ".rds"))
}

FillNaEmpty <- function(element){
    if(length(element)==0) element <- NA
    return(element)
}

MakeSpecificCrisprsetList <- function(res.file.df, res.file.path, condition.label){
  crisprset.list <- apply(res.file.df, MARGIN = 1, function(res.row){
    if(res.row["sample.label"] %in% c(condition.label) & file.exists(file.path(res.row["aligner.res.dir.path"], res.file.path))){
      detected.variant.info <- readRDS(file.path(res.row["aligner.res.dir.path"], res.file.path))
      variant.n <- length(detected.variant.info) / 4 # variant number
      mut.label.vec <- character()
      mut.cnt.vec <- numeric()
      mut.consensus.vec <- character()
      mut.insertion.vec <- character()
      mut.deletion.vec <- character()
      mut.mh.vec <- character()
      for(variant.ind in 0:(variant.n - 1)){
          # TODO:ここからコンセンサス配列とカウント情報を抜き出す、倍数を取得する必要あり
          mut.lebel <- detected.variant.info[variant.ind * 4 + 1]$mut.name
          if(length(grep("S", mut.lebel)) > 0) next # ignore substitution
          if(length(grep("M", mut.lebel)) > 0) next # ignore large indel

          mut.cnt.vec <- c(mut.cnt.vec, detected.variant.info[variant.ind * 4 + 2]$mut.cnt)
          mut.insertion <- as.character(detected.variant.info[variant.ind * 4 + 4]$indel.seq$insert.seq)
          mut.insertion.vec <- c(mut.insertion.vec, FillNaEmpty(mut.insertion))
          mut.deletion <- as.character(detected.variant.info[variant.ind * 4 + 4]$indel.seq$deletion.seq)
          mut.deletion.vec <- c(mut.deletion.vec, FillNaEmpty(mut.deletion))
          mut.mh <- detected.variant.info[variant.ind * 4 + 4]$indel.seq$mh.seq.vec
          if(length(mut.mh) == 0) mut.mh <- ""
          mut.mh.vec <- c(mut.mh.vec, mut.mh)
          
          # modify label
          if(length(grep("I", mut.lebel)) > 0){
              mut.lebel <- paste0(mut.lebel ,":", mut.insertion) # add inseretion seq
              mut.lebel <- gsub("^[-0-9]+:", ":", mut.lebel)
          }
          if(nchar(mut.mh) > 0){
              mut.lebel <- paste0(mut.lebel, ":mh", nchar(mut.mh)) # add mh deletion length
          }
          mut.label.vec <- c(mut.label.vec, mut.lebel)
      }
      variant.df <- data.frame(mut.label=mut.label.vec, mut.cnt=mut.cnt.vec, mut.insertion=mut.insertion.vec, mut.deletion=mut.deletion.vec, mut.mh=mut.mh.vec)
      res <- list(variant.df)
      names(res) <- res.row["sample.name"]
    }else if(res.row["sample.label"] %in% c(condition.label)){
      res <- NULL
    }else{
      res <- NULL
    }
    return(res)
  })
  return(Filter(Negate(is.null), crisprset.list)) # return list without NULL element
}

CalcKLValue <- function(x.cigar_freqs, y.cigar_freqs){
  # prepare vector
  variants.names.vec <- unique(c(rownames(x.cigar_freqs), rownames(y.cigar_freqs)))
  # browser()

  # NOTE:
  # add small pseudocounts of 0.5 in one sample but not the other for avoinding division by zero. This method was inspired by [Felicity Allen, et.al, 2017]

  x.variants.vec <- rep(0.0, length(variants.names.vec))
  names(x.variants.vec) <- variants.names.vec
  x.variants.vec[rownames(x.cigar_freqs)] <- x.variants.vec[rownames(x.cigar_freqs)] + as.numeric(x.cigar_freqs[,1])
  for(cigar.name in names(x.variants.vec)){
    if(x.variants.vec[cigar.name] == 0){
        x.variants.vec[cigar.name] = x.variants.vec[cigar.name] + 0.5
    }
  }
  x.variants.vec <- x.variants.vec / sum(x.variants.vec)

  y.variants.vec <- rep(0.0, length(variants.names.vec))
  names(y.variants.vec) <- variants.names.vec
  y.variants.vec[rownames(y.cigar_freqs)] <- y.variants.vec[rownames(y.cigar_freqs)] + as.numeric(y.cigar_freqs[,1])
  for(cigar.name in names(y.variants.vec)){
    if(y.variants.vec[cigar.name] == 0){
        y.variants.vec[cigar.name] = y.variants.vec[cigar.name] + 0.5
    }
  }
  y.variants.vec <- y.variants.vec / sum(y.variants.vec)

  return(sum(x.variants.vec * log(x.variants.vec / y.variants.vec)))
}

# detected pattern vs predicted pattern
CalcKLdivergence <- function(x.cigar_freqs, y.indel.data, can_classify_by_range=FALSE, x.label, y.label, filepath){
  if(is.null(x.cigar_freqs) & is.null(y.indel.data)){
    return(0)
  }else if(is.null(x.cigar_freqs) | is.null(y.indel.data)){
    return(Inf)
  }
  AggregateVariantsFreq <- function(cigar_freqs){
    aggregated.freq.table <- aggregate(reads ~ label, data = data.frame(label = gsub("^[-0-9]+:", ":", rownames(cigar_freqs)),  reads = as.numeric(cigar_freqs[,1])), sum)
    position.independent.cigar_freqs <- matrix(
      aggregated.freq.table$reads,
      ncol = 1
    )
    rownames(position.independent.cigar_freqs) <- aggregated.freq.table$label
    return(position.independent.cigar_freqs)
  }
  x.indel.data <-  AggregateVariantsFreq(x.cigar_freqs)
  y.indel.data <- AggregateVariantsFreq(data.frame(y.indel.data[,1]))

  if(x.label==y.label){
      SaveTable(x.indel.data, paste0(filepath, "_CRISPResso2_indel_", x.label))
      SaveTable(y.indel.data, paste0(filepath, "_InDelPhi_indel_", y.label))
  }

  if(can_classify_by_range){
      x.indel.data <- matrix(position.independent.x.cigar_freqs[which(rownames(position.independent.x.cigar_freqs) %in% rownames(y.indel.data)),], ncol = 1)
      rownames(x.indel.data) <- rownames(position.independent.x.cigar_freqs)[which(rownames(position.independent.x.cigar_freqs) %in% rownames(y.indel.data))]
      x.indel.class.vec = c(
        x.indel.data["1I",1],
        x.indel.data["1D",1],
        x.indel.data["2D",1],
        sum(x.indel.data[,1]) - (x.indel.data["1I",1] + x.indel.data["1D",1] + x.indel.data["2D",1])
      )
      names(x.indel.class.vec) <- c("1I", "1D", "2D", ">2D")
      x.indel.data <- DataFrame(x.indel.class.vec)

      y.indel.class.vec = c(
        y.indel.data["1I",1],
        y.indel.data["1D",1],
        y.indel.data["2D",1],
        100 - (y.indel.data["1I",1] + y.indel.data["1D",1] + y.indel.data["2D",1])
      )
      names(y.indel.class.vec) <- c("1I", "1D", "2D", ">2D")
      y.indel.data <- DataFrame(y.indel.class.vec)
  }
  return(
    CalcKLValue(x.indel.data, y.indel.data) +
    CalcKLValue(y.indel.data, x.indel.data)
  )
}

# predicted indel pattern vs detected indel pattern
MakeHeatMap <- function(fltr.res.file.df, predicted.data.list, y.axis.title, x.axis.title, file.name.excluding.ext, condition.label, class_type){
  
  res.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅰ]crispresso_alignment_data/[Ac]crispresso_nhej/variants.count.total.list.rds"
  detected.data.list <- MakeSpecificCrisprsetList(fltr.res.file.df, res.path, condition.label)

  # Replace name
  names(predicted.data.list) <- gsub("-[A-Za-z0-9]*$", paste0("-", condition.label), names(predicted.data.list))

  # make matrix for plot
  
  predicted.comp.kl.mat <- matrix(NA, ncol = length(predicted.data.list))
  colnames(predicted.comp.kl.mat) <- names(predicted.data.list)
  rownames(predicted.comp.kl.mat) <- "temp"
  for(ind in 1:length(detected.data.list)){
    predicted.comp.kl.vec <- numeric(0)
    machiato.nhej.crisprset.name <- names(detected.data.list[[ind]])
    for(temp.name in names(predicted.data.list)){
      detected.data.df <- getElement(detected.data.list[[ind]], machiato.nhej.crisprset.name)
      
      # extrast only single deletions and one-base insertion data
      detected.data.unpredict.range.index <- !grepl(":", detected.data.df$mut.label)
      # TODO:ここから一処理ずつ訂正していく
      for(insertion_size in 2:200){
          detected.data.unpredict.range.index <- detected.data.unpredict.range.index | grepl(paste0(":", insertion_size, "I"), detected.data.df$mut.label) # 1bp> insertions
      }
      for(deletion_size in 61:200){
          detected.data.unpredict.range.index <- detected.data.unpredict.range.index | grepl(paste0(":", deletion_size, "D"), detected.data.df$mut.label) # 60bp> deletions
      }
      detected.data.unpredict.range.index <- detected.data.unpredict.range.index | grepl(paste0(","), detected.data.df$mut.label) # multiple indels
      # remove substitution and complex vatiant
      detected.data.predict.range.df <- DataFrame(machiato.nhej=detected.data.df[!detected.data.unpredict.range.index,]$mut.cnt)
      rownames(detected.data.predict.range.df) <- detected.data.df[!detected.data.unpredict.range.index,]$mut.label
     # Calculate KL divergence
     if(class_type == "base"){
         predicted.comp.kl.vec <- c(predicted.comp.kl.vec, CalcKLdivergence(detected.data.predict.range.df, getElement(predicted.data.list, temp.name), can_classify_by_range=FALSE, names(detected.data.list[[ind]]), temp.name, file.name.excluding.ext))
     }else if(class_type == "range"){
         predicted.comp.kl.vec <- c(predicted.comp.kl.vec, CalcKLdivergence(detected.data.predict.range.df, getElement(predicted.data.list, temp.name), can_classify_by_range=TRUE, names(detected.data.list[[ind]]), temp.name, file.name.excluding.ext))
     }

      if(length(predicted.comp.kl.vec) < 2){
        names(predicted.comp.kl.vec) <- c(temp.name)
      }else{
        names(predicted.comp.kl.vec) <- c(names(predicted.comp.kl.vec)[1:length(predicted.comp.kl.vec)-1], temp.name)
      }
    }
    predicted.comp.kl.mat <- rbind(predicted.comp.kl.mat, predicted.comp.kl.vec)
    # just below sentence was written due to misunderstanding relationship between colnames and rownames.
    # actually, I should change colnames because colnames mean predicted score. but I changed rownames.
    # I need a few change to fix it. but it works correctly without dealing with it. So it remains.
    rownames(predicted.comp.kl.mat) <- c(rownames(predicted.comp.kl.mat)[1:nrow(predicted.comp.kl.mat)-1], gsub(paste0("-", condition.label, "$"), "", c(machiato.nhej.crisprset.name)))
  }
  # rownames(predicted.comp.kl.mat) # MaChIAto
  # colnames(predicted.comp.kl.mat) # Predicted Score
  predicted.comp.kl.mat <- predicted.comp.kl.mat[-1, ]
  predicted.comp.kl.mat <- predicted.comp.kl.mat[, which(colnames(predicted.comp.kl.mat) %in% paste0(rownames(predicted.comp.kl.mat), "-", condition.label))]

  reorder.ind.vec <- unlist(sapply(paste0(rownames(predicted.comp.kl.mat), "-", condition.label), function(x){
    return(which(colnames(predicted.comp.kl.mat) %in% x))
  }))
  predicted.comp.kl.mat <- predicted.comp.kl.mat[, reorder.ind.vec]
  
  predicted.comp.kl.p.mat <- predicted.comp.kl.mat # just renamed same matrix
  colnames(predicted.comp.kl.p.mat) <- gsub(paste0("-", condition.label, "$"), "", colnames(predicted.comp.kl.p.mat))
  
  # make heat map
  options(warn=-1)
  # browser()
  indelphi.comp.kl.p <- ggplot(melt(log2(predicted.comp.kl.p.mat[,ncol(predicted.comp.kl.p.mat):1]))
    , aes(x=Var1, y=Var2, fill=value)) +
  geom_raster(aes(fill = value)) +
  labs(title ="Symmetrized KL divergence"
    , x = x.axis.title
    , y = y.axis.title) +
  theme(axis.ticks = element_blank()
    , axis.text.x = element_text(angle = 270, hjust = 0, size=7)
    , axis.text.y = element_text(size=7)
    , plot.title = element_text(size = 10)
    , axis.title.x = element_text(size = 10)
    , axis.title.y = element_text(size = 10)
  ) +
  scale_fill_distiller(palette = "Spectral", name = "Symmetrized KL Divergence (log2)", limits = c(min(log2(predicted.comp.kl.p.mat)), max(log2(predicted.comp.kl.p.mat))))
  # scale_fill_distiller(palette = "Spectral", name = "Symmetrized KL Divergence (log2)", limits = c(min(log2(predicted.comp.kl.p.mat)), max(log2(predicted.comp.kl.p.mat))))

  SaveTable(predicted.comp.kl.p.mat, file.name.excluding.ext)
  ggsave(file = paste0(file.name.excluding.ext, ".png")
    , plot = indelphi.comp.kl.p, dpi = 350, width = 11, height = 10)


  # make boxplot
  match.value.vec <- numeric(0)
  nonmatch.value.vec <- numeric(0)
  for(row.ind in 1:nrow(predicted.comp.kl.mat)){
    for(col.ind in 1:ncol(predicted.comp.kl.mat)){
      if(row.ind == col.ind){
        match.value.vec <- c(match.value.vec, log2(predicted.comp.kl.mat)[row.ind, col.ind])
      }else{
        nonmatch.value.vec <- c(nonmatch.value.vec, log2(predicted.comp.kl.mat)[row.ind, col.ind])
      }
    }
  }

  boxplot.df <- data.frame(variable = c(rep("corresponding.target", length(match.value.vec)), rep("noncorresponding.target", length(nonmatch.value.vec))),
    value = c(match.value.vec, nonmatch.value.vec))
  
  options(warn = -1)
  ylim.vec <- c(min(boxplot.df$value), max(boxplot.df$value))
  box.plot <- ggplot(boxplot.df, aes(x=variable, y=value, fill = factor(variable))) + 
    geom_boxplot(width = 0.5, alpha = 0.8, show.legend = FALSE, outlier.shape = NA) + 
    geom_jitter(width = 0.05, alpha=0.4, colour = "#F29F05") +
    coord_cartesian(ylim = ylim.vec) +
    scale_fill_manual(values = c("#0528F2", "#F24405")) +
    theme(axis.ticks = element_blank()
      , axis.text.x = element_text(size=15)
      , axis.text.y = element_text(size=15)
      , axis.title.y = element_text(size = 20)
      , legend.position = "none"
    ) +
    labs(x = "", y = "Symmetrized KL Divergence (log2)") +
    scale_y_continuous(breaks=seq(round(ylim.vec[1] - 1), round(ylim.vec[2] + 1), 1))
  
  ggsave(file = paste0(file.name.excluding.ext, "_boxplot.png"), plot = box.plot, dpi = 350, width = 10, height = 10)
  options(warn = 0)

  # run test
  res.t.test <- t.test(match.value.vec, nonmatch.value.vec, paired = FALSE)
  sink(paste(file.name.excluding.ext, "_ttest.txt"))
  print("t-test")
  if(!is.null(res.t.test)){
    print(paste0("t-statistic : ", res.t.test$statistic))
    print(paste0("degrees of freedom : ", res.t.test$parameter))
    print(paste0("p-value : ", res.t.test$p.value))
  }else{
    print("t-test was failed...")
  }
  print(paste0("number of corresponding.target : ", length(match.value.vec)))
  print(paste0("mean of corresponding.target : ", mean(match.value.vec)))
  print(paste0("SD of corresponding.target : ", sd(match.value.vec)))
  print(paste0("number of noncorresponding.target : ", length(nonmatch.value.vec)))
  print(paste0("mean of noncorresponding.target : ", mean(nonmatch.value.vec)))
  print(paste0("SD of noncorresponding.target : ", sd(nonmatch.value.vec)))
  sink()

  options(warn=0)
}

plot.diagonal.heatmap <- function(data, len, min.value, max.value, file.name.excluding.ext, gene.name, title.name){
  diagonal.kl.p <- ggplot(melt(data[,len:1]), aes(x=Var1, y=Var2, fill=value)) +
    geom_raster(aes(fill = value)) +
    labs(title = title.name , x = "Condition", y = gene.name) +
    theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 270, hjust = 0, size=7), axis.text.y = element_text(size=7), plot.title = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) +
    scale_fill_distiller(palette = "Spectral", name = "Symmetrized KL Divergence", limits = c(min.value, max.value))
  SaveTable(data, file.name.excluding.ext)
  ggsave(file = paste0(file.name.excluding.ext, ".png")
    , plot = diagonal.kl.p, dpi = 350, width = 16, height = 2)
}

getAverExp <- function(sc, gene=rownames(sc@data), conds, rest=F){
  mtExp <- as.matrix(expm1(sc@data[gene, ]))
  avgExp <- list()
  avgExp[[conds]] <- unlist(rowMeans(mtExp[, WhichCells(sc, ident = conds)]))
  if(rest)
    avgExp[["rest"]] <- unlist(rowMeans(mtExp[, WhichCells(sc, ident.remove = conds)]))
  return(avgExp)
}

getPlot <- function(sc, rep, col, key, scale=F){
  df <- as.data.frame(sc@dr[[rep]]@cell.embeddings[, 1:2])
  df[, "key"] <- key
  df[, "val"] <- col
  my.x <- colnames(df)[1]
  my.y <- colnames(df)[2]
  
  p <- ggplot(data=df, aes(x=df[,my.x], y=df[,my.y], key=df[, "key"])) +
    theme(line = element_blank(),
          title = element_blank(),
          axis.text = element_blank(),
          axis.ticks =  element_blank(),
          axis.line =  element_blank())
  #         text = element_blank(),
  #         title = element_blank())
  
  if(scale)
    p <- p + scale_colour_gradient2(low = ("blue"), mid = "yellow", high = ("red"), midpoint = median(df$val), name="Legend")
  
  return(p)
}

runFindAllMks <- function(sc){
  result <- list()
  df <- FindAllMarkers(sc)
  df <- df %>% group_by(cluster) %>% arrange(desc(avg_logFC), .by_group = T)
  avrExp <- lapply(levels(sc@ident), function(x) getAverExp(sc, gene = df$gene, conds = x, rest = T))
  names(avrExp) <- levels(sc@ident)
  cluster_exp <- NULL
  rest_exp <- NULL
  for(i in levels(df$cluster)){
    my.df <- df %>% filter(cluster==i)
    cluster_exp <- c(cluster_exp, avrExp[[i]][[1]][my.df$gene])
    rest_exp <- c(rest_exp, avrExp[[i]][[2]][my.df$gene])
  }
  df$cluster_exp <- cluster_exp
  df$rest_exp <- rest_exp
  df <- df %>% group_by(cluster) %>% arrange(desc(avg_logFC), .by_group = T)
  df <- as.data.frame(df)
  result[["all"]] <- df
  df <- df[, c(5,2,7,6,8,9)]
  df <- df %>% filter(p_val_adj<0.05)
  df <- as.data.frame(df)
  result[["filter"]] <- df
  return(result)
  
}


runFindMks <- function(sc, ident1, ident2 = NULL){
  result <- list()
  df <- FindMarkers(sc, ident.1 = ident1, ident.2 = ident2)
  df$gene <- rownames(df)
  df <- df %>% arrange(desc(avg_logFC))
  if(!is.null(ident2)){
    avrExp <- lapply(c(ident1, ident2), function(x) getAverExp(sc, gene = df$gene, conds = x))
    names(avrExp) <- c(ident1, ident2)
    cluster_exp <- list()
    for(i in names(avrExp)){
      cluster_exp[[i]] <- avrExp[[i]][[1]][df$gene]
    }
    df$cluster_exp <- cluster_exp[[1]]
    df$rest_exp <- cluster_exp[[2]]
  }else{
    avrExp <- getAverExp(sc, gene = df$gene, conds = ident1, rest = T)
    df$cluster_exp <- avrExp[[1]][df$gene]
    df$rest_exp <- avrExp[[2]][df$gene]
  }
  df <- df %>% arrange(desc(avg_logFC))
  result[["all"]] <- df
  df <- df[, c(5,2,6,7,8)]
  df <- df %>% filter(p_val_adj<0.05)
  df <- as.data.frame(df)
  result[["filter"]] <- df
  return(result)
}

color.list <- c(brewer.pal(12, "Paired"),brewer.pal(12, "Set3"),brewer.pal(8, "Pastel2"),colorRampPalette(c("grey20","grey70"))(4))
color.gradient <- function(x, colors=c("blue","yellow","red"), colsteps=50, dscale=NULL) {
  if (!is.null(dscale)) {
    x[which(x < dscale[1])] <- dscale[1]
    x[which(x > dscale[2])] <- dscale[2]
  }
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x,na.rm = T),max(x,na.rm = T), length.out=colsteps)) ] )
}

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "enrichResult"),
          function(x, showCategory = 30, color = "p.adjust", layout = "kk", ...) {
              emapplot.enrichResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "gseaResult"),
          function(x, showCategory = 30, color = "p.adjust", layout = "kk", ...) {
              emapplot.enrichResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })

##' @rdname emapplot
##' @exportMethod emapplot
setMethod("emapplot", signature(x = "compareClusterResult"),
          function(x, showCategory = 30, color = "p.adjust", layout = "kk", ...) {
              emapplot.compareClusterResult(x, showCategory = showCategory,
                                    color = color, layout = layout, ...)
          })
		  
		  
		  
##' @rdname emapplot
##' @importFrom igraph graph.empty
##' @importFrom igraph add_vertices
##' @importFrom igraph graph.data.frame
##' @importFrom igraph delete.edges
##' @importFrom igraph V "V<-"
##' @importFrom igraph E "E<-"
##' @importFrom reshape2 melt
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_gradientn
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_void
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @author Guangchuang Yu
emapplot.enrichResult <- function(x, showCategory = 30, color="p.adjust", layout = "kk", ...) {
    n <- update_n(x, showCategory)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    if (is.numeric(n)) {
        y <- y[1:n,]
    } else {
        y <- y[match(n, y$Description),]
        n <- length(n)
    }


    if (n == 0) {
        stop("no enriched term found...")
    } else if (n == 1) {
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- y$Description
        V(g)$color <- "red"
        return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
    } else {
        id <- y[,1]
        geneSets <- geneSets[id]

        n <- nrow(y) #
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y$Description

        for (i in 1:n) {
            for (j in i:n) {
                w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }

        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
        g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
        g <- delete.edges(g, E(g)[wd[,3] < 0.2])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))

        cnt <- sapply(geneSets[idx], length)
        V(g)$size <- cnt

        colVar <- y[idx, color]
        V(g)$color <- colVar
    }


    p <- ggraph(g, layout=layout)

    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    }

    p + geom_node_point(aes_(color=~color, size=~size)) +
        geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
        scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
        ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
        scale_size(range=c(3, 8))
}

geneID.compareClusterResult <- function(x) {
    as.character(x@compareClusterResult$geneID)
}
geneInCategory.compareClusterResult <- function(x) {
    setNames(strsplit(geneID(x), "/", fixed = TRUE), x@compareClusterResult$ID)

}
merge_compareClusterResult <- function(yy) {
	yy_union<- yy[!duplicated(yy$ID),]
	rownames(yy_union) <- yy_union$ID
	#做一个for循环来合并基因
	for(i in seq_len(dim(yy)[1])) {
	    gene1 <- unlist(strsplit(yy[i,"geneID"],"/", fixed = TRUE))
		gene2 <- unlist(strsplit(yy_union[yy[i,"ID"],"geneID"],"/", fixed = TRUE))
		genes <- union(gene1,gene2)
		gene_ID_uni <- paste(genes,collapse="/")
		yy_union[yy[i,"ID"],"geneID"] <- gene_ID_uni	
	}
	return(yy_union)
}

emapplot.compareClusterResult <- function(x, showCategory = 30, color="p.adjust", layout = "kk", ...) {

n <- enrichplot:::update_n(x, showCategory)
geneSets <- geneInCategory(x) 
y <- as.data.frame(x)
  if (is.numeric(n)) {
        y <- y[1:n,]
    } else {
        y <- y[match(n, y$Description),]
        n <- length(n)
    }

#geneSets <- geneInCategory(x) ## use core gene for gsea result

#进行数据结构转换，将同个ID（Description）的基因进行合并


y_union <- merge_compareClusterResult(y)

    if (n <= 1) {
        stop("no enriched term found...")
    } else {
	    id <- y_union$ID
		geneSets <- geneSets[id]
        n <- nrow(y_union) #
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y_union$Description

        for (i in 1:n) {
            for (j in i:n) {
                w[i,j] = enrichplot:::overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }

        wd <- melt(w)
        wd <- wd[wd[,1] != wd[,2],]
        wd <- wd[!is.na(wd[,3]),]
		g <- graph.data.frame(wd[,-3], directed=FALSE)
        E(g)$width=sqrt(wd[,3] * 5)
		#先不对边的阈值进行筛选
        g <- delete.edges(g, E(g)[wd[,3] < 0.05])
        idx <- unlist(sapply(V(g)$name, function(x) which(x == y_union$Description)))

        cnt <- sapply(geneSets[idx], length)
        V(g)$size <- cnt

        colVar <- y_union[idx, color]
        V(g)$color <- colVar
	}
	

#下面把饼图添上去就行了	
#先画出所有需要的饼图，使用的数据仅为两列
data_pie <- as.matrix(y[,1:2])
ID_unique <- unique(data_pie[,2])
Cluster_unique <- unique(data_pie[,1])
ID_Cluster_mat <- matrix(0,length(ID_unique),length(Cluster_unique))
rownames(ID_Cluster_mat) <- ID_unique
colnames(ID_Cluster_mat) <- Cluster_unique
for(i in seq_len(dim(data_pie)[1])) {
    ID_Cluster_mat[data_pie[i,2],data_pie[i,1]] <- ID_Cluster_mat[data_pie[i,2],data_pie[i,1]]+1
}

dff <- melt(t(ID_Cluster_mat))
names(dff) <- c("Cluster","ID","value")


#先画出每一张饼图
nn <- dim(ID_Cluster_mat)[2]
pies <- rep(1,nrow(dff)/nn)
k <-  1
for(i in seq(1,nrow(dff),nn)) {
     df <- dff[i:(i+nn-1),]
	 outfile <- paste0(df[1,"ID"], ".png")
	 ggplot(df, aes(x = 1, value, fill = Cluster)) +
        geom_col() + coord_polar(theta = 'y') +
        ggtitle(df[1,"ID"]) +
        theme_void() + theme_transparent() +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(6), hjust = 0.5))+
			  ggsave(outfile, bg = "transparent")
	pies[k] <- outfile
	k <- k + 1
 
}

#将饼图加到网络图中
	  p <- ggraph(g, layout=layout)

    if (length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
    }


	
	
	
aa <- p$data

radius <- sqrt(aa$size / 3000)
	

p + geom_image(aes(image = pies,x=aa$x,y=aa$y,size = I(radius)))+
scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
theme_void()
}

emapplot.compareClusterResult(xx)


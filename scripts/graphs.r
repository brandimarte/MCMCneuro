library(gRbase)

myiplot <- function (g, ...)
{
  V(g)$size <- 13
  V(g)$label.cex <- 2
  plot (g,...)
}

gr = as.matrix(read.table ("adjM6HPp1Met3.dat", check.names=FALSE))
gG1 <- as(gr, "igraph")

# Ajeite o grafo do jeito que preferir e MANTENHA a janela aberta
tkplot(gG1)
xy <- tkplot.getcoords(1)
png ("g6_p1.png", width = 900, height = 900)
myiplot (gG1, layout=xy)
dev.off ()

gr = as.matrix(read.table ("adjM6HPp3Met3.dat", check.names=FALSE))
gG2 <- as(gr, "igraph")
png ("g6_p3.png", width = 900, height = 900)
myiplot (gG2, layout = xy)
dev.off ()

# --------------------------------------------
# Short instructions for installing Rgraphviz:

# Install Graphviz version >= 2.2
# source("http://bioconductor.org/biocLite.R")
# biocLite("RBGL")
# biocLite("graph")
# biocLite("Rgraphviz")

# Citations to these packages:
# citation("graph")
# citation("RBGL")
# citation("Rgraphviz")


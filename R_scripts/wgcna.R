library(WGCNA)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(tcltk)

adj <- TOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj)
network <- simplify(network)  # removes self-loops
V(network)$color <- blocks$colors
V(network)$name <- names(blocks$colors)
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
plot(network, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)

group_ids = unique(blocks$colors)
group_color <- brewer.pal(length(group_ids), 'Set1')
group_color_fill <- paste0(group_color, '20')

plot(network, vertex.color = blocks$colors, layout=layout.fruchterman.reingold(network), 
     vertex.size = 10, edge.arrow.size = 0.2,
     edge.color = rgb(0.5, 0.5, 0.5, 0.2),
     mark.col = group_color_fill,
     mark.border = group_color)

tkplot(network)

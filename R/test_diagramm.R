nodes <-
  create_node_df(
    n = 6,
    type = c("a","b","b","c","d","d"),
    label= c("insect sting reaction in the ER",
             "anaphylaxis",
             "local reaction",
             "treatment",
             "acute serum samples\nwithin 2h of reaction",
             "baseline serum samples\nafter > 2 weeks"
             ),

    #color = c("red", "green",
    #          "grey", "blue"),
    #value = c(3.5, 2.6, 9.4, 2.7),
    shape = "rectangle",
    width = c(2,1,1,1,2,2),
    x = c(0, -.6, .6, -0.6,-0, -0),
    y = c(0,-0.5,-0.5,  -1,-1.5,-2),
    height = c(0.3,0.3,0.3,0.3,0.4,0.4))


edges <- create_edge_df(from = c(1,1,2,4,5,3),
                        to =   c(2,3,4,5,6,5),
                        headport = c(5,0,0,0,0,0))
# Add the node data frame to the
# graph object to create a graph
# with nodes
gr <- create_graph() %>%
  add_node_df(node_df = nodes) %>%
  add_edge_df(edge_df = edges) %>%
  add_global_graph_attrs(
    attr = "ortho",
    value = 1,
    attr_type = "graph") %>%
  add_global_graph_attrs(value="black",
                         attr = "fontcolor",
                         attr_type = "node") #%>% render_graph()
  export_graph(file_name = "analysis/figures/flow.png", file_type = "png",
               title = NULL,
               width = NULL,
               height = NULL)

  library(network)
  library(sna)
  n <- network(rgraph(6, tprob = 0.2), directed = FALSE)

library(ggnetwork)

  df <- ggnetwork(n, layout = "fruchtermanreingold", cell.jitter = 0.75)
  df$x <- c(0, -.6, .6, -0.6,-0, -0)
  ggplot(df,
         aes(x = x, y = y, xend = xend, yend = yend))+
   #ggplot(n,aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(color = "grey50") +
    theme_blank()

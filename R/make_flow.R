### mAKE THE DIAGRAMM FOR THE fLOW CHART

make_flow <- function(){
require(DiagrammeRsvg)
require(DiagrammeR)

# Create a node data frame
nodes <-
  create_node_df(
    n = 7,
    type = "a",
    label= c(paste0("All cases in the\nEuropean Anaphylaxis Registry\n",length(data4$b_submitdate)),
             paste0("Formal anaphylaxis definition met\n",
                    length(which(data4$reaction_type_brown=="anaphylaxis"))),
             paste0("Insects as elicitors\n",
                    length(which(data4$d_elicitor_gr5== "insects" &
                                   data4$reaction_type_brown=="anaphylaxis"))),
             paste0("Other elicitors\n",
                    length(which(data4$d_elicitor_gr5!= "insects"&
                                   data4$reaction_type_brown=="anaphylaxis"))),
             paste(paste0("Yellow-jackets = ",
                          length(which(data4$q_340_insects== "yellow jacket"&
                                         data4$reaction_type_brown=="anaphylaxis"))),
                   paste0("Bees = ",
                          length(which(data4$q_340_insects== "bee"&
                                         data4$reaction_type_brown=="anaphylaxis"))),
                   paste0("Hornets = ",
                          length(which(data4$q_340_insects== "hornet"&
                                         data4$reaction_type_brown=="anaphylaxis"))),
                   paste0("Other insects = ",
                          length(which(data4$d_insect_gr4=="other" &
                                         data4$reaction_type_brown=="anaphylaxis"))),
                   sep="\n"
             ),
             paste0("Elicitor unknown\n",
                    length(which(data4$d_elicitor_gr5=="unkown" &
                                   data4$reaction_type_brown=="anaphylaxis"))),
             paste0("Other known\nelicitors\n",
                    length(which(data4$d_elicitor_gr5!="insects"&
                                   data4$d_elicitor_gr5!="unkown" &
                                   data4$reaction_type_brown=="anaphylaxis")))
    ),

    #color = c("red", "green",
    #          "grey", "blue"),
    #value = c(3.5, 2.6, 9.4, 2.7),
    shape= "rectangle",
    width = c(2.5,3,1.8,1.2,1.8,1.2,1.2),
    x = c(2, 2, 1.1, 3.1,1.1,4.7,3.1),
    y = c(0,-1,-2,-2,-3,-2,-3),
    height = c(rep(0.6,4),1,0.6,0.6))

edges <- create_edge_df(from = c(1,2,2,3,4,4),
                        to =   c(2,3,4,5,6,7))
# Add the node data frame to the
# graph object to create a graph
# with nodes
create_graph() %>%
  add_node_df(node_df = nodes) %>%
  add_edge_df(edge_df = edges) %>%
  add_global_graph_attrs(
    attr = "splines",
    value = 1,
    attr_type = "graph") %>%
  add_global_graph_attrs(value="black",
                         attr = "fontcolor",
                         attr_type = "node") %>% #render_graph()
  export_graph(file_name = "analysis/figures/flow.png", file_type = "png",
               title = NULL,
               width = 1300,
               height = 1000)
}

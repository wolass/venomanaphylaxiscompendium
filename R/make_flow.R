### mAKE THE DIAGRAMM FOR THE fLOW CHART

make_flow <- function(){
require(DiagrammeRsvg)
require(DiagrammeR)

# Create a node data frame
nodes <-
  create_node_df(
    n = 8,
    type = "a",
    label= c(paste0("All cases in the\nEuropean Anaphylaxis Registry\n",
                    format(length(data4$b_submitdate),
                           big.mark = ",")),
             paste0("Formal anaphylaxis definition met\n",
                    format(length(which(data4$reaction_type_brown=="anaphylaxis")),
                           big.mark = ",")),
             paste0("Sex- and age-matching of venom induced\nanaphylaxis and other known elicitors\n",
                    format(nrow(age_sex_matched %>%
                                  filter(grouping == "other")),
                           big.mark = ",")),
             paste0("Insects as elicitors\n\n",
                    format(length(which(data4$d_elicitor_gr5== "insects" &
                                   data4$reaction_type_brown=="anaphylaxis")),
                           big.mark = ",")),
             paste0("Other elicitors\n(matched control group)\n",
                    format(nrow(age_sex_matched %>%
                                  filter(grouping == "other")),
                           big.mark = ",")),
                    #format(length(which(data4$d_elicitor_gr5!= "insects"&
                    #               data4$reaction_type_brown=="anaphylaxis")),
                    #       big.mark = ",")),
             paste(paste0("yellow-jackets = ",
                          format(length(which(data4$q_340_insects== "yellow jacket"&
                                         data4$reaction_type_brown=="anaphylaxis")),
                                 big.mark = ",")),
                   paste0("bees = ",
                          format(length(which(data4$q_340_insects== "bee"&
                                         data4$reaction_type_brown=="anaphylaxis")),
                                 big.mark = ",")),
                   paste0("hornets = ",
                          length(which(data4$q_340_insects== "hornet"&
                                         data4$reaction_type_brown=="anaphylaxis"))),
                   paste0("other insects = ",
                          length(which(data4$d_insect_gr4=="other" &
                                         data4$reaction_type_brown=="anaphylaxis"))),
                   sep="\n"
             ),
             paste0("Elicitor\nunknown\n",
                    format(length(which(data4$d_elicitor_gr5=="unkown" &
                                   data4$reaction_type_brown=="anaphylaxis")),
                           big.mark = ",")),
             paste0("food = ",
                    format(age_sex_matched %>%
                                  count(d_elicitor_gr5) %>%
                             pull(n),
                           big.mark = ",")[1],
                    "\ndrugs = ",
                    format(age_sex_matched %>%
                             count(d_elicitor_gr5) %>%
                             pull(n),
                           big.mark = ",")[2],
                    "\nother = ",
                    format(age_sex_matched %>%
                             count(d_elicitor_gr5) %>%
                             pull(n),
                           big.mark = ",")[4]
             )

    ),

    #color = c("red", "green",
    #          "grey", "blue"),
    #value = c(3.5, 2.6, 9.4, 2.7),
    shape= "rectangle",
    x =      c(   2,     2,   2, 1.1, 3.1, 1.1,    4.2, 3.1),
    y =      c(-0.5, -1.25,  -2,-2.9,-2.9,  -4,  -1.25,  -4),
    width =  c(   3,     3,   3, 1.8, 1.8, 1.8,      1, 1.8),
    height = c(  .6,    .6,  .6,  .8,  .8, 0.8,    0.6, 0.8))

edges <- create_edge_df(from = c(1,2,3,3,4,2,5),
                        to =   c(2,3,4,5,6,7,8))
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
  add_global_graph_attrs(value="white",
                         attr = "fillcolor",
                         attr_type = "node") %>% #render_graph()
  export_graph(file_name = "analysis/figures/flow.png", file_type = "png",
               title = NULL,
               width = 1000,
               height = 900)
}

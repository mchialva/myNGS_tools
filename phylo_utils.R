# A wrapper around this nice script to collect group coordinates to plot phylogenetic tree annotation using geom_segment() and ggtree()
# Credits: https://rawgit.com/valentinitnelav/valentinitnelav.github.io/master/assets/2018-01-07-ggtree/2018-01-07-ggtree.html#51_annotate_groups_with_segments_and_text
# Usage coord_gr(tree, offset_y)
# tree: ggtree-linke phylogenetic tree
# offset_y: a value indicating the reverse of the distance between categories

coord_gr<-function(tree, offset_y)
{
  tree_dt<-data.table(tree$data)[isTip == TRUE][order(y)]
  coord_groups <- tree_dt[, .(y1 = y[1],
                            y2 = y[.N],
                            angle = mean(angle),
                            n = .N), # optional - helps with counting
                        by = .(group, 
                               id_gr = rleid(group, 
                                             prefix = "grp"))]
coord_groups[, y_mid := rowMeans(.SD), .SDcols = c("y1", "y2")]
coord_groups[, y1_adj := ifelse(y1 == y2, y1 +offset_y, y1-offset_y)]
coord_groups[, y2_adj := ifelse(y1 == y2, y2 -offset_y, y2+offset_y)]
coord_groups[, angle_adj := ifelse(angle %between% c(90, 180), 
                                   yes = angle + 180,
                                   no = ifelse(angle > 180 & angle <= 270,
                                               yes = angle - 180,
                                               no = angle))]

coord_groups[, hjust_adj := ifelse(angle %between% c(90, 270), yes = 1L, no = 0L)]
return(list("tree_dt"=tree_dt, "coord_groups"=coord_groups))
}

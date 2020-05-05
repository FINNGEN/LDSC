#example script for using the ComplexHeatmap package to plot correlations

library(ComplexHeatmap)
library(circlize)
library(jpeg)

  #prep input for plotting
  ##input_dataset: first column is list of UKB phenotypes, second column is list of UKB categories, remaining columns are correlations with each FG drug
prepare_dataset_for_heatmap <- function(input_dataset){
  rownames(input_dataset) <- input_dataset[,1]
  input_dataset[,c(1,2)] <- NULL
  input_dataset <- t(input_dataset)
}

fg_ukb_rg_matrix_pthres_subset_final_for_heatmap <- prepare_dataset_for_heatmap(fg_ukb_rg_matrix_pthres_subset_final)

  #set column split as as column with UKB phenotype categories
UKB_split=fg_ukb_rg_matrix_pthres_subset_final$CATEGORY 
  #set row split as column with drug categories
drugs_split=drugs$category 
  #set main heatmap legend parameters
heatmap_legend_param=list(col_fun = col_fun, 
                          title = "Correlation", 
                          title_gp = gpar(fontsize=12, fontface = "bold"),
                          at = c(-1, -0.5, 0, 0.5, 1), 
                          legend_height = unit(10, "cm")) 
  #set row_title for all rows
row_title_general=c("FinnGen Drugs")
  #set column_title for all columns
column_title_general=c("UKB Phenotypes")
  #set row_title for row categories
row_title="%s"
  #set column_title for column categories
column_title="%s"
  #set color scheme
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

jpeg("heatmap_clustering_by_UKB_and_drug_final_subset.jpeg", width = 2400*2, height = 2700, type="cairo")
heatmap = Heatmap(fg_ukb_rg_matrix_pthres_subset_final_for_heatmap, 
             col=col_fun,
             width = unit(120, "cm"), height = unit(40, "cm"),
             
             cluster_rows=T, clustering_distance_rows='euclidean', clustering_method_rows='complete',
             show_row_dend=F, 
             cluster_columns=T, clustering_distance_columns = 'euclidean', clustering_method_columns = 'complete',
             show_column_dend = F,
             
             column_names_gp=gpar(fontsize=10),
             column_names_side="bottom",
             column_names_max_height=unit(60, "cm"),
             column_names_rot=90,
             
             row_names_gp=gpar(fontsize=15),
             row_names_side="left",
             row_names_max_width=unit(100, "cm"), 
            
             heatmap_legend_param=heatmap_legend_param,
             
             column_split=UKB_split,
             column_gap = unit(4, "mm"),
             column_title_rot = 90,
             column_title=column_title,
             column_title_side="top",
             column_title_gp=gpar(fontface="bold",fontsize=15),
             
             row_split=drugs_split,
             row_gap = unit(4, "mm"),
             row_title_rot=0,
             row_title=row_title,
             row_title_side="right",
             row_title_gp=gpar(fontface="bold",fontsize=15),
             
             border=T)
draw(heatmap,
     
     row_title=row_title_general, 
     row_title_side="left",
     row_title_gp=gpar(fontface="bold",fontsize=20),
     
     column_title=column_title_general, 
     column_title_side="bottom", 
     column_title_gp=gpar(fontface="bold",fontsize=20),
     
     use_raster=T, raster_device="jpeg", raster_quality=3)
dev.off()

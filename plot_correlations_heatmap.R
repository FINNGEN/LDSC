#import rg, pvals
rg_pvals <- rg_pvals.txt

#import phenotypes 
phenos <- phenolist.txt

#cat rg, pvals with phenotypes
merge <- cbind(phenos, rg_pvals)

#get drug names
ATC_translate_minimal <- read_delim("~/ATC_translate_minimal.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
ATC_translate_minimal_subset <- subset(ATC_translate_minimal, select=c("Id", "Longname"))

#merge with drug names
merge_final <- merge(merge, ATC_translate_minimal_subset, by.x="name_of_FG_drug_col", by.y="Id") 

#create and save heatmap 
rg_heatmap <- ggplot(data = merge_final, aes(x=FG_drug, y=UKB_pheno, fill=genetic_correlation)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) + 
  geom_tile() +  coord_fixed() + theme_minimal() +
  geom_text(aes(label = p_val), color = "black", size = 2.5) + theme(axis.text.x = element_text(angle = 45,
                                                                                                size = 10, vjust=1, hjust=1)) +
  labs(x = "FinnGen Drug\n", y = "UKBiobank Phenotype\n") 
ggsave(plot = rg_heatmap, height = 10, width = 10*3, dpi = 300, filename = "rg_heatmap.jpeg", path="~/rg")

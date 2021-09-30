library(tidyverse)
library(ggpubr)
library(rstatix)
library(cowplot)
library(picante)
library(phyloseq)
library(microbiome)
library(vegan)
library(microbiomeutilities)
library(ANCOMBC)
library(ape)
library(phangorn)
library(GUniFrac)


# filter dataframe by rownames
filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
}

#########################    IMPORTING DATA #########################  
split_taxonomy <- function(otu) {
  # select taxa column and split taxonomy column into separate columns using ; delimiter
  taxonomy_phylo <- dplyr::select(otu, "taxonomy") %>%
    separate("taxonomy", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")
  return(taxonomy_phylo)
}

read_tab_delim <- function(df) {
  # read all tab delimited files using these params
  df_out <- read.table(df, row.names = 1, header = 1, sep="\t", comment.char = "", check.names = F )
  return(df_out)
}

load_phylo <- function(otu, taxa, mapping, tree=NULL) {
  # convert to phyloseq and return list 
  phylo_otu <- otu_table(otu, taxa_are_rows = T)
  
  phylo_tax <- tax_table(as.matrix(taxa))
  
  phylo_map <- sample_data(mapping)
  
  if (exists("tree")){
    phylo_tree <- read_tree(tree)
    return(merge_phyloseq(phylo_otu, phylo_tax, phylo_tree, phylo_map))
  }
  else {
    return(merge_phyloseq(phylo_otu, phylo_tax, phylo_map))
  }
}

prep_4_phyloseq <- function(otu, mapping, tree=NULL){
  # read files
  otu_taxa <- read_tab_delim(otu)
  metadata <- read_tab_delim(mapping)
  # convert otu into matrix and drop taxonomy col
  count_matrix <- as.matrix(subset(otu_taxa, select=-taxonomy))
  # split taxonomy
  taxonomy_phylo <- split_taxonomy(otu_taxa)
  # make sure taxa table row names match otu 
  rownames(taxonomy_phylo) <- rownames(otu_taxa)
  # use load phylo function to convert all files to phyloseq objs
  if (exists("tree")){
    # if user wants to input a phylogenetic tree this code is called
    out <- load_phylo(count_matrix, taxonomy_phylo, metadata, tree = tree)
  }
  else {
    # without a tree
    out <- load_phylo(count_matrix, taxonomy_phylo, metadata)
  }
  return(out)
}

######################### NORMALISATION #########################

# normalise data with minimum sum scaling 
transform_mss <- function(x) {
  x_untransformed <- x
  if (length(is(x)) == 1 && class(x) == "phyloseq"){
    x <- abundances(x)
  }
  else {
    print("not a phyloseq or data not integers, exiting")
    stop()
  }
  xt <- t(min(colSums(x)) * t(x) / colSums(x))
  otu_table(x_untransformed)@.Data <- xt
  x_out <- x_untransformed
  return(x_out)
}

######################### PLOTTING AND STATS #########################

# calculate adonis R2 and p-value
phyloseq_adonis <- function(ps, dist_matrix, factor, ...) {
  
  
  metadata_df <- as(sample_data(ps), "data.frame") 
  
  if(sum(is.na(metadata_df[[factor]])) > 0){
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {
    ps_ad <- adonis(dist_matrix ~ metadata_df[[factor]],
                    data = metadata_df, ...)
    return(ps_ad)
  }
}


# create statistical formula 
xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

# plot boxplot with stats
plot_boxplot <- function(df, variable_col, value_col, 
                         fill_var="fill", comparisons_list, xlab, ylab, 
                         p_title=NULL, multiple_groups=FALSE, col_palette,
                         group.order=NULL, ...){
  
  if (!is.null(group.order)){
    print(group.order)
    df[, variable_col] <- factor(df[, variable_col], levels = group.order)
    print(df)
  }
  
  formula <- xyform(value_col, variable_col)
  
  if (multiple_groups == TRUE) {
    stat_variance <- df %>% 
      kruskal_test(formula)
    stat_test <- df %>%  
      pairwise_wilcox_test(formula, comparisons = comparisons_list,
                           p.adjust.method = "BH") %>% 
      add_significance() %>% 
      add_xy_position(x = variable_col) %>% 
      filter(p.adj < 0.05)
  } 
  else if (multiple_groups == FALSE) {
    stat_test <- df %>% 
      wilcox_test(formula) %>% 
      add_significance() %>% 
      add_xy_position(x = variable_col) %>% 
      filter(p < 0.05)
  }
  
  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(df, aes_string(x=variable_col, y=value_col, 
                                fill = variable_col, color=variable_col )) + 
    geom_boxplot(color="black", alpha=0.8,outlier.shape=5, outlier.size=1) + 
    geom_point(size = 1.5, position = position_jitterdodge()) + 
    labs(x=xlab, y=ylab) + 
    stat_boxplot(color="black", geom = "errorbar", width=0.2) 
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot + 
    theme_cowplot() + 
    ggtitle(p_title) + 
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(), 
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = col_palette) + 
    scale_color_manual(values = col_palette) + 
    rotate_x_text(angle = 45)
  
  if (dim(stat_test)[1] == 0){
    plot_out <- final_plot
  }
  else {
    plot_out <- final_plot +
      stat_pvalue_manual(stat_test, label = "p.adj.signif", hide.ns = T,
                         inherit.aes = FALSE, ...)
  }
  
  return(plot_out)
}

# plot scatter plot with correlation if desired

#TODO - rewrite so it can be used for betadiv as well 
plot_scatter <- function(df, x, y, point_color, 
                         line_color, fill_color, 
                         xlabel, ylabel, 
                         corr.method=NULL, ...) {
  p <- ggplot(data=df, mapping = aes(x = .data[[ x ]], y = .data[[ y ]] )) +
    geom_point(aes(color=point_color), size=2.5) +
    geom_smooth(method = "lm", color=line_color, fill=fill_color) +
    theme_bw() + 
    theme(legend.position = "None", 
          axis.title.x = element_text(size=14),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=12)) + 
    xlab(xlabel) + 
    ylab(ylabel)
  
  if (!is.null(corr.method)) {
    p <- p + stat_cor(method=corr.method, ...)
    return(p)
  }
  else {
    return(p)
  }
}

######################### ALPHA DIVERSITY #########################

# calculate shannon effective
Shannon.E <- function(x) {
  summed <- sum(x)
  shannon.e <- round(exp(-sum(x[x>0]/summed*log(x[x>0]/summed))), digits = 2)
  return(shannon.e)
}

# calculate faiths phylogenetic diversity
Faiths <- function(x, tree) {
  pd <- pd(x, tree, include.root = F)
  return(pd[,1])
}

# calculate richness
Richness <- function(x) {
  observed <- sum(x[x>0.5]^0)
  return(observed)
}

# calculate all alpha diversity matrices and return dataframe 
calc_alpha <- function(ps) {
  mat_in <- abundances(ps) %>% 
    t()
  
  tree <- phy_tree(ps)
  
  diversity <- setNames(data.frame(matrix(ncol = 3, nrow = nsamples(ps))), 
                        c("Richness", "Shannon.Effective", "Faiths.PD"))
  rownames(diversity) <- rownames(meta(ps))
  
  diversity$Richness <- apply(mat_in, 1, Richness)
  diversity$Shannon.Effective <- apply(mat_in, 1, Shannon.E)
  diversity$Faiths.PD <- Faiths(mat_in, tree)
  
  return(diversity)
}



######################### BETA DIVERSITY #########################
### wrap these into two functions 

# Calculate Gunifrac in phyloseq objects 
phyloseq_gunifrac <- function(phyloseq_obj, asdist=TRUE) {
  # input phyloseq obj
  # returns dist matrix of GUnifrac distance which can 
  # be used as input for ordinate() function 
  
  #extract otu table and tree
  otu_tab <- t(otu_table(phyloseq_obj))
  tree <- phy_tree(phyloseq_obj)
  
  # root tree at midpoint
  rooted_tree <- midpoint(tree)
  
  # calculate gunifrac
  gunifrac <- GUniFrac(otu_tab, rooted_tree, 
                       alpha = c(0.0,0.5,1.0))$unifracs
  # alpha param is weight on abundant lineages, 0.5 has best power
  # so we will extract this 
  if (asdist == T){
    return(stats::as.dist(gunifrac[, , "d_0.5"]))
  }
  else{
    return(gunifrac[, , "d_0.5"])
  }
}


calc_betadiv <- function(ps, ord_method, dist) {
  if (ord_method %in% c("NMDS", "MDS", "PCoA")) {
    
    
    if (distance %in% c("unifrac", "wunifrac", "bray")) {
      
      dist_metric <- distance(ps, dist)
      ord <- ordinate(ps, ord_method, dist_metric)
      
      return_list <- list("Distance_Matrix" = dist_metric, 
                          "Ordination" = ord$points)
      return(return_list)
    }
    else if (distance %in% c("gunifrac")) {
      
      gu <- phyloseq_gunifrac(ps)
      ord <- ordinate(ps, ord_method, gu)
      
      return_list <- list("Distance_Matrix" = gu, 
                          "Ordination" = ord$points)
      return(return_list)
    }
  else {
    print("Ordination method not supported, supported methods are:
          NMDS, MDS, PCoA")
    stop()
  }
} 
}

plot_beta_div <- function(ps, ordination, dist_matrix, group_variable, add_ellipse=FALSE) {
  # add support for colors and elipses
  # significance
  ad <- phyloseq_adonis(ps, dist_matrix, group_variable)
  
  
  plot <- plot_ordination(ps, ordination , color=group_variable)
  plot$layers[[1]] <- NULL
  
  plot_out <- plot + geom_point(size=3, alpha=0.75) +
    theme_bw() + 
    scale_fill_manual(values=cols) +
    scale_color_manual(values = cols) +
    labs(caption = bquote(Adonis~R^2 ~ .(round(ad$aov.tab$R2[1], 2)) ~~ p-value ~ .(ad$aov.tab$`Pr(>F)`[1])))
  
  return(plot_out)
}



######################### Differential Abundance #########################


pairwise_wilcox_asv <- function(asv_tab, factor, case, control, 
                                sample_col, meta, p_adjust="BY"){
  
  # initialise pval matrix 
  p.val <- matrix(NA, nrow=nrow(asv_tab), ncol=1, 
                  dimnames=list(row.names(asv_tab), "p.val"))
  
  featmat <- as.matrix(asv_tab)
  
  for (row in rownames(featmat)){
    x <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==control) %>% pull(.data[[ sample_col ]])]
    y <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==case) %>% pull(.data[[ sample_col ]])]
    
    
    p.val[row, ] <- wilcox.test(x, y, exact=FALSE)$p.value
  }
  
  p.val.adj <- p.val %>% 
    as.data.frame() %>% 
    adjust_pvalue("p.val", method = p_adjust)
  
  return(p.val.adj)
}

ancom_da <- ancombc(phyloseq = ps, formula = "Group", 
                    p_adj_method = "BH", zero_cut = 0.33, 
                    group = "Group", struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

ancom_res_df <- data.frame(
  Species = row.names(ancom_da$res$beta),
  beta = unlist(ancom_da$res$beta),
  se = unlist(ancom_da$res$se),
  W = unlist(ancom_da$res$W),
  p_val = unlist(ancom_da$res$p_val),
  q_val = unlist(ancom_da$res$q_val),
  diff_abn = unlist(ancom_da$res$diff_abn))

fdr_ancom <- ancom_res_df %>%
  dplyr::filter(q_val < 0.05)

dim(fdr_ancom)

da_ancom <- ancom_res_df %>% 
  filter(diff_abn == T) 
da_asvs <- da_ancom$Species

taxonomy <- as.data.frame(tax_table(ps))
da_taxonomy <- taxonomy %>%
  select(best_hit)%>% 
  filter(best_hit %in% da_asvs)

da_ancom <- da_ancom %>% 
  rownames_to_column(var="Group") %>% 
  column_to_rownames("Species") %>% 
  merge(da_taxonomy, by=0)

da_ancom <- da_ancom %>% 
  mutate(enriched_in = ifelse(beta > 0, "Stroke", 
                              "Sham"))
# order by FC
da_ancom_sort <- da_ancom %>% 
  arrange(beta) %>% 
  mutate(best_hit=factor(best_hit, levels=best_hit))


#Differential abundance 


ggplot(da_ancom_sort, aes(x = best_hit, y = beta, color = enriched_in)) +
  geom_point(size = 5) +
  labs(y = "\nLog2 Fold-Change - Stroke vs. Sham", x = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "none") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dotted") + 
  scale_color_manual(values=c("#122C34", "#ED254E"))
ggsave("figures/DA_taxa_LFC.pdf", device = "pdf", dpi=300, height = 10, width = 14)
```

import seaborn as sns 
import pandas as pd 
import numpy as np
import zipfile 
from skbio.stats import distance
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None

## FUNCTIONS

# calculate unpaired wilcoxon - only required if manually adding p-vals 
def wilcox_test(df, metric, group_col, groups):
    x = df[df[group_col] == groups[0]][metric]
    y = df[df[group_col] == groups[1]][metric]
    return mannwhitneyu(x, y)  

# get identifier from unzipped qiime2 artifact
def get_zip_name(zipobj):
    # get path from zipfile
    fullpath = zipobj.namelist()[1]
    return(fullpath.split("/")[0])

def plot_alpha(alphadiv, alpha_metrics, group, cols, labels, pairs):
    for i in alpha_metrics:
        # boxplot
        p = sns.boxplot(x=alphadiv[group], y=alphadiv[i], width=0.4, palette=cols)
        p.set_ylabel(labels[i])
        # calc and add stats
        annotator = Annotator(p, pairs, data=alphadiv, x=group, y=i)
        annotator.configure(test='Mann-Whitney', text_format='star')
        annotator.apply_and_annotate()
        # save fig to pdf
        plt.savefig(i + ".pdf", dpi=300, bbox_inches="tight")
        # show figure inline
        plt.show()

def ordinate(dissimilarity_matrix, metadata, group, ord_method="NMDS",random_seed=42):
    if ord_method in ["NMDS", "PCoA", "MDS"]:
        if ord_method == "NMDS":
            # initiate MDS object
            comp = MDS(n_components=2, dissimilarity="precomputed", metric=False, random_state=random_seed)
            # calculate ordination 
            ordination = pd.DataFrame(comp.fit_transform(dissimilarity_matrix), columns=["NMDS1", "NMDS2"], index=dissimilarity_matrix.index)
        elif ord_method in ["PCoA", "MDS"]:
            # initiate MDS object
            comp = MDS(n_components=2, dissimilarity="precomputed", metric=True, random_state=random_seed)
            # calculate ordination 
            ordination = pd.DataFrame(comp.fit_transform(dissimilarity_matrix), columns=["MDS1", "MDS2"], index=dissimilarity_matrix.index)
    else:
        print("Ordination method not supported")
        return     
    # join with metadata for plotting
    metadata_plot = pd.concat([ordination, metadata], axis=1)
    
    
    # reorder metadata as permanova requires dataframes in same order 
    metadata_ord = metadata.reindex(dissimilarity_matrix.index)
    # create distance matrix class from gunifrac matrix
    dm = distance.DistanceMatrix(dissimilarity_matrix,ids=dissimilarity_matrix.index)
    # calculate permanova
    p_res = distance.permanova(dm, metadata_ord, column=group)
    # currently permanova implementation in scikit-bio does not calculate r2 so we do this manually
    r2 = round(1 - 1 / (1 + p_res[4] * p_res[3] / (p_res[2] - p_res[3] - 1)), 3)
    # append r2 value to permanova results
    p_res["R2"] = r2
    # return metadata combined with ordination results and permanova output
    return metadata_plot, p_res

def plot_beta_div(plot_df, permanova_res, group, cols, ord_method="NMDS"):
    # define axes names
    if ord_method == "NMDS":
        axes = ["NMDS1", "NMDS2"]
    elif ord_method in ["PCoA", "MDS"]:
        axes = ["MDS1", "MDS2"]
      
    # save annotation text with pval and r2 to simplify code below
    annotation =  "r2= " + str(permanova_res["R2"]) + "\n" + "p= "  + str(permanova_res["p-value"])
    # plot NMDS and add p-value and r2
    p = sns.scatterplot(data=plot_df, x=axes[0], y=axes[1],hue="Group", palette=cols, s=80, alpha=0.8)
    p.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # add label with x,y coordinates
    x_coord = max(plot_df[axes[0]]) * 0.75
    y_coord = min(plot_df[axes[1]]) * 1.75
    # add annotations and save/show figure
    plt.text(x_coord, y_coord, annotation , horizontalalignment='left', size='small', color='black', style="italic", weight="regular")
    plt.savefig("ordination_plot.pdf", dpi=300, bbox_inches="tight")
    plt.show()
    

    
def format_taxonomy(taxonomy, col):
    # split into columns for each rank
    taxonomy[["Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"]] = taxonomy[col].str.split(';',expand=True)
	# replace empty cells with nan
    taxonomy.replace("", np.nan, inplace=True)
	# get empty indices
    taxonomy_empty = taxonomy.isnull()
	# forward fill na
    taxonomy_format =  taxonomy.fillna(method='ffill', axis=1)
	# using indices from before add prefix to ffilled cells
    taxonomy_format.update('unknown_' + taxonomy_format.mask(~taxonomy_empty))
	# combine formatted taxonomy into one new column, overwriting old taxonomy
    taxonomy_format["taxonomy"] = taxonomy_format[taxonomy_format.columns[1:]].apply(lambda x: ";".join(x.dropna().astype(str)),
    axis=1)
	# drop extra cols and return
    taxonomy_format.drop(["Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species", col], axis=1, inplace=True)

    return taxonomy_format

def add_taxonomy(ancom_res, taxonomy):
    ancom_res_taxa = pd.merge(ancom_res, taxonomy, left_index=True, right_index=True)
    return ancom_res_taxa


# this code needs improved
def highest_classified(ancom_res, taxcol):
    # get and add highest classified taxonomy 
    tax_format = format_taxonomy(ancom_res[[taxcol]], col=taxcol)
    # get genus level taxonomy
    tax_format["highest_classified"] = tax_format["taxonomy"].str.split(';',expand=True)[5]
    ancom_res["highest_classified"] = tax_format["highest_classified"]
    return ancom_res

def pseudo_asv(ancom_res):
    pseudo_asvs = [f'ASV_{i}' for i in range(1, len(ancom_res.index) + 1)]
    ancom_res.index = pseudo_asvs
    return ancom_res
    
def plot_ancom(ancom, metadata, taxonomy,groups, cols, taxcol="Taxon"):
    # add taxonomy 
    ancom_tax = add_taxonomy(ancom, taxonomy)
    # convert q-val and FC to numeric (changed after merging with taxonomy)
    ancom_tax.iloc[:,[0,4]] =  ancom_tax.iloc[:,[0,4]].apply(pd.to_numeric)
    # get sig 
    ancom_sig = ancom_tax[ancom_tax["q-value"] < 0.05]
    
    ancom_sig_tax = highest_classified(ancom_sig, taxcol)
    
    # convert beta (log FC) to log2 FC
    ancom_sig_tax["Log2FC"] = np.log2(np.exp(ancom_sig_tax["beta"]))
    ancom_sig_tax["Enriched"] = np.where(ancom_sig_tax["Log2FC"] < 0, groups[0], groups[1])
    
    # sort by Log2FC and make names unique 
    ancom_sig_tax = pseudo_asv(ancom_sig_tax)
    ancom_sig_tax_sort = ancom_sig_tax.sort_values(by ='Log2FC' , ascending=False)
    
    ancom_sig_tax_sort["ASV"] = ancom_sig_tax_sort.index + ";" + ancom_sig_tax_sort["highest_classified"]
    p = sns.stripplot(data=ancom_sig_tax_sort, x="Log2FC", y="ASV", hue="Enriched", palette=cols, s=10)
    p.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    p.set_xlabel("Log2 Fold-Change " + str(groups[0]) + " vs " + str(groups[1]))
    plt.savefig("ancom_da_plot.pdf", dpi=300, bbox_inches="tight")
    plt.show()
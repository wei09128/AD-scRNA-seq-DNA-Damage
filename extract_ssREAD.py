import pandas as pd
import glob
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# 1. R Logic:
r_logic = """
library(qs)
library(Seurat)

process_and_extract <- function(file_path, target_genes) {
    obj <- qs::qread(file_path)
    
    if ("RNA" %in% names(obj@assays)) {
        obj@assays$RNA@meta.features <- data.frame(row.names = rownames(obj@assays$RNA))
    }
    
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
    
    n_pcs <- min(20, ncol(obj) - 1)
    obj <- RunUMAP(obj, dims = 1:n_pcs, verbose = FALSE)
    
    umap_df <- as.data.frame(Embeddings(obj, reduction = "umap"))
    
    existing_genes <- intersect(target_genes, rownames(obj))
    if (length(existing_genes) > 0) {
        gene_expr <- FetchData(obj, vars = existing_genes)
    } else {
        gene_expr <- data.frame(NoMatch = rep(0, ncol(obj)), row.names = colnames(obj))
    }
    
    meta <- obj@meta.data
    meta$source_file <- basename(file_path)
    
    return(list(umap = umap_df, expr = gene_expr, meta = meta))
}
"""

robjects.r(r_logic)
process_func = robjects.globalenv['process_and_extract']

h_genes_list = ["ATM", "TP53", "CHEK2", "PARP1", "CDKN1A", "CDKN2A", "RELA", "BRCA1", "ABL1"]
m_genes_list = [g.capitalize() for g in h_genes_list] 
all_genes_search = robjects.StrVector(h_genes_list + m_genes_list)

def safe_convert(r_obj):
    if isinstance(r_obj, pd.DataFrame): return r_obj
    with localconverter(robjects.default_converter + pandas2ri.converter):
        return pandas2ri.rpy2py(r_obj)

files = glob.glob("*.qsave")

for f in files:
    base_name = os.path.splitext(f)[0] # Get filename without .qsave
    print(f"--- Processing: {f} ---")
    
    try:
        results_r = process_func(f, all_genes_search)
        
        umap_pd = safe_convert(results_r[0])
        expr_pd = safe_convert(results_r[1])
        meta_pd = safe_convert(results_r[2])
        
        combined_temp = pd.concat([umap_pd, expr_pd, meta_pd], axis=1)
        
        found_h = [g for g in h_genes_list if g in combined_temp.columns]
        found_m = [g for g in m_genes_list if g in combined_temp.columns]

        # Species Logic
        combined_temp['h_sum'] = combined_temp[found_h].sum(axis=1) if found_h else 0
        combined_temp['m_sum'] = combined_temp[found_m].sum(axis=1) if found_m else 0
        h_count = (combined_temp[found_h] > 0).sum(axis=1) if found_h else 0
        m_count = (combined_temp[found_m] > 0).sum(axis=1) if found_m else 0

        is_human = (combined_temp['h_sum'] > combined_temp['m_sum']) | \
                   ((combined_temp['h_sum'] == combined_temp['m_sum']) & (h_count >= m_count) & (h_count > 0))
        
        is_mouse = (combined_temp['m_sum'] > combined_temp['h_sum']) | \
                   ((combined_temp['m_sum'] == combined_temp['h_sum']) & (m_count > h_count))

        # 1. Save Mouse CSV for this specific file
        mouse_cells = combined_temp[is_mouse].copy()
        if not mouse_cells.empty:
            mouse_cells = mouse_cells.drop(columns=found_h + ['h_sum', 'm_sum'], errors='ignore')
            mouse_cells.to_csv(f"{base_name}_mouse.csv")
            print(f"   >> Saved {base_name}_mouse.csv ({len(mouse_cells)} cells)")

        # 2. Save Human CSV for this specific file
        human_cells = combined_temp[is_human].copy()
        if not human_cells.empty:
            human_cells = human_cells.drop(columns=found_m + ['h_sum', 'm_sum'], errors='ignore')
            human_cells.to_csv(f"{base_name}_human.csv")
            print(f"   >> Saved {base_name}_human.csv ({len(human_cells)} cells)")

    except Exception as e:
        print(f"  !! Error processing {f}: {e}")

print("\nAll qsave files have been converted to individual CSVs.")
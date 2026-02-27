import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import os

# Define the genes for both species
GENE_DICT = {
    "human": ["ATM", "TP53", "CHEK2", "PARP1", "CDKN1A", "CDKN2A", "RELA", "BRCA1", "ABL1"],
    "mouse": ["Atm", "Trp53", "Chek2", "Parp1", "Cdkn1a", "Cdkn2a", "Rela", "Brca1", "Abl1"]
}

def find_col(df, possible_names):
    for name in possible_names:
        found = next((c for c in df.columns if c.lower() == name.lower()), None)
        if found: return found
    return None

# Loop through all files in the current directory
for filename in os.listdir('.'):
    if filename.endswith(".csv"):
        print(f"--- Processing: {filename} ---")
        
        # Load Data once per CSV
        df = pd.read_csv(filename, index_col=0, low_memory=False)
        base_name = os.path.splitext(filename)[0]

        # Identify UMAP and Cluster columns
        u1 = find_col(df, ['UMAP_1', 'umap_1', 'rna.umap_1'])
        u2 = find_col(df, ['UMAP_2', 'umap_2', 'rna.umap_2'])
        cluster_col = find_col(df, ['predicted.id', 'seurat_clusters', 'cell_type'])

        if not u1 or not u2:
            print(f"Skipping {filename}: UMAP coordinates not found.")
            continue

        # 1. TRY BOTH SPECIES FOR GENE EXPRESSION
        for species, genes_to_show in GENE_DICT.items():
            valid_genes = [g for g in genes_to_show if g in df.columns]

            if not valid_genes:
                print(f"   - No {species} genes found in {filename}.")
                continue

            print(f"   - Generating {species} plot...")
            n_genes = len(valid_genes)
            cols = 4
            rows = math.ceil(n_genes / cols)
            fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4))
            
            # Handle single-gene edge case
            if n_genes == 1: axes = [axes]
            else: axes = axes.flatten()

            for i, gene in enumerate(valid_genes):
                plot_df = df.sort_values(by=gene)
                sc = axes[i].scatter(plot_df[u1], plot_df[u2], c=plot_df[gene], cmap='YlOrRd', s=1.5, alpha=0.8)
                axes[i].set_facecolor('#f0f0f0')
                axes[i].set_title(f"{gene}", fontsize=14, fontweight='bold')
                axes[i].axis('off')
                plt.colorbar(sc, ax=axes[i], fraction=0.046, pad=0.04)

            # Clean up empty subplots
            for j in range(i + 1, len(axes)): axes[j].axis('off')

            plt.suptitle(f"{species.upper()} DNA Damage Signaling: {base_name}", fontsize=22, y=0.98)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            
            # SAVE NAME: filename_human.png or filename_mouse.png
            plt.savefig(f"{base_name}_{species}.png", dpi=300)
            plt.close()

        # 2. GENERATE CELL TYPE PLOT (One per CSV)
        import seaborn as sns
        import matplotlib.pyplot as plt
        
        # 1. Define a high-contrast palette manually for maximum distinction
        # This uses a mix of bold primary and secondary colors
        high_contrast_colors = [
            "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", 
            "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", 
            "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", 
            "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"
        ]
        
        # Inside your loop, replace the cell type plot section with this:
        if cluster_col:
            print(f"   - Plotting cell types with high-contrast colors...")
            plt.figure(figsize=(14, 9))
            
            # Get unique categories and sort them for a consistent legend
            categories = sorted(df[cluster_col].unique())
            
            # Apply the high-contrast palette
            sns.scatterplot(
                data=df, 
                x=u1, y=u2, 
                hue=cluster_col, 
                hue_order=categories,
                s=3,                # Increased size slightly for better color pop
                palette=high_contrast_colors[:len(categories)], 
                alpha=1.0,          # Full opacity to distinguish colors better
                edgecolor=None
            )
            
            plt.title(f"Cell Types: {base_name}", fontsize=18, fontweight='bold')
            
            # Enhanced Legend
            plt.legend(
                bbox_to_anchor=(1.02, 1), 
                loc='upper left', 
                markerscale=3,      # Large legend dots
                frameon=True,       # Add a border to the legend box
                fontsize=11
            )
            
            plt.savefig(f"{base_name}_celltypes.png", dpi=300, bbox_inches='tight')
            plt.close()

print("\nAll files processed successfully.")
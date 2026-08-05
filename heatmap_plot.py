# This script reads the rbp_family_matrix.csv and mainly plots the heatmap using matplotlib and seaborn

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Load data
df = pd.read_csv('rbp_family_matrix.csv', index_col=0)

# Custom sorting to make the graph more readable
TOP_N_RBPS = 40        
TOP_M_FAMILIES = 15    

# Calculate totals to identify the most-studied targets across the literature
top_rbps = df.sum(axis=1).nlargest(TOP_N_RBPS).index
top_families = df.sum(axis=0).nlargest(TOP_M_FAMILIES).index

# Subset the matrix to high-yield intersections before log transformation
df_filtered = df.loc[top_rbps, top_families]
log_matrix = np.log1p(df_filtered)

# Define custom color palette
colors = ["#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac", "#053061"]
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_virus", colors)

sns.set(font_scale=0.8)

# Plotting with sns
g = sns.clustermap(
    log_matrix.T,
    cmap=custom_cmap,
    figsize=(18, 8),  # 22 > 18
    xticklabels=True,
    yticklabels=True,
    dendrogram_ratio=0.001,  # Removed it later
    cbar_pos=(0.02, 0.75, 0.01, 0.15),
    cbar_kws={"label": "log(1 + paper count)"}
)

# Remove dendrograms
g.ax_row_dendrogram.set_visible(False)
g.ax_col_dendrogram.set_visible(False)
g.cax.set_visible(False)


g.ax_heatmap.set_position([0.10, 0.22, 0.74, 0.72])
g.ax_heatmap.yaxis.tick_left()
g.ax_heatmap.yaxis.set_label_position("left")
g.ax_heatmap.tick_params(axis='y', labelleft=True, labelright=False, pad=2, length=0)
g.ax_heatmap.set_xlabel('RBP', fontname='Arial', fontsize=14, labelpad=10)
g.ax_heatmap.set_ylabel('\nVirus families', fontname='Arial', fontsize=14, labelpad=10)

plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontname='Arial', fontsize=14)
plt.setp(g.ax_heatmap.get_yticklabels(), fontstyle='italic', fontname='Arial', fontsize=12)


# Extract the reordered row and column indices calculated by the dendrogram
reordered_row_indices = g.dendrogram_row.reordered_ind
reordered_col_indices = g.dendrogram_col.reordered_ind

clustered_output_matrix = log_matrix.T.iloc[reordered_row_indices, reordered_col_indices]



log_matrix.T.to_csv('rbp_family_log_transformed_filtered.csv')
clustered_output_matrix.to_csv('rbp_family_matrix_CLUSTERED_ORDER_filtered.csv')




plt.show()
plt.close()
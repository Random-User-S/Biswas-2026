# This script plots the upset plot

# pip install openpyxl for reading the xlsx


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# File paths
excel_path = "rbp_master_literature_matrix.xlsx"
csv_func_path = "rbp_functional_classes_matrix.csv"

# Metadata & classifications
# Animal viruses only
virus_families = [
    'Adenoviridae', 'Arenaviridae', 'Bunyavirales', 'Caliciviridae', 'Coronaviridae',
    'Deltavirus', 'Filoviridae', 'Flaviviridae', 'Hepadnaviridae', 'Hepeviridae',
    'Herpesviridae', 'Orthomyxoviridae', 'Papillomaviridae', 'Paramyxoviridae',
    'Parvoviridae', 'Picornaviridae', 'Polyomaviridae', 'Poxviridae', 'Reoviridae',
    'Retroviridae', 'Rhabdoviridae', 'Togaviridae'
]

func_classes = [
    'RNA splicing', 'Translation', 'RNA stability/decay', 'RNA processing',
    'RNA transport', 'Gene regulation', 'Stress granules', 'Innate immunity',
    'RNA modification', 'RNA helicase'
]

# Color palette for functional classes
func_palette = {
    "RNA splicing":        "#7FB3D5",
    "RNA processing":      "#E67E7E",
    "RNA transport":       "#B48ECA",
    "RNA stability/decay": "#8FD0C8",
    "RNA modification":    "#D7CF4F",
    "RNA helicase":        "#55B8C7",
    "Translation":         "#A8D08D",
    "Gene regulation":     "#D8A27A",
    "Stress granules":     "#D7B6D8",
    "Innate immunity":     "#A8A58A",
}
# (https://imagecolorpicker.com/)


left_set_bar_color = "#2F6FAF"

# Load datasets
df_lit = pd.read_excel(excel_path)
df_func = pd.read_csv(csv_func_path)


lit_rbp_col = "Query_RBP" if "Query_RBP" in df_lit.columns else "RBP"
func_rbp_col = "RBP" if "RBP" in df_func.columns else "RBP_Symbol"

# Merge datasets
df_merged = pd.merge(df_lit, df_func, left_on=lit_rbp_col, right_on=func_rbp_col, how="inner")

# Standardize RBP column name to prevent KeyError in downstream grouping
df_merged['RBP'] = df_merged[lit_rbp_col]

# Convert indicators to binary integers
df_bin_virus = (df_merged[virus_families] > 0).astype(int)
df_bin_func = df_merged[func_classes].isin(['+', '✓', 1, True]).astype(int)

df_merged['Num_Active_Viruses'] = df_bin_virus.sum(axis=1)
df_merged['Num_Func_Classes'] = df_bin_func.sum(axis=1)

tot_papers_col = "Total_Papers_All_Families" if "Total_Papers_All_Families" in df_merged.columns else "Total Papers (All Families)"
df_merged['Total_Papers'] = df_merged[tot_papers_col]

# Filter virus families (keep families with at least 5 RBPs; arbitrary value)
set_sizes_series = df_bin_virus.sum(axis=0).sort_values(ascending=False)
set_sizes_series = set_sizes_series[set_sizes_series >= 5]

active_v_sets = set_sizes_series.index.tolist()
total_v_sizes = [int(df_bin_virus[v].sum()) for v in active_v_sets]

# Group RBPs by virus family intersection combination
df_merged['v_combo'] = df_bin_virus.apply(lambda r: tuple(r[v] for v in active_v_sets), axis=1)

combo_df = df_merged.groupby('v_combo').agg(
    Intersection_Size=('RBP', 'count'),
    Num_Active_Viruses=('v_combo', lambda x: sum(x.iloc[0])),
    Num_Func_Classes=('Num_Func_Classes', 'sum'),
    Max_Papers=('Total_Papers', 'max')
).reset_index()

combo_df = combo_df[combo_df['v_combo'].apply(lambda c: any(c))]



# Sorting logic:
# 1. Has functional classes (No_Functions == 0 first)
# 2. Number of active virus families (descending)
# 3. Total functional classes (descending)
# 4. Max papers published (descending)
combo_df['No_Functions'] = (combo_df['Num_Func_Classes'] == 0).astype(int)

sorted_combo_df = combo_df.sort_values(
    by=['No_Functions', 'Num_Active_Viruses', 'Num_Func_Classes', 'Max_Papers'],
    ascending=[True, False, False, False]
).reset_index(drop=True)

top_combo_df = sorted_combo_df.head(25).copy()
num_intersections = len(top_combo_df)
x_coordinates = np.arange(num_intersections)

all_func_breakdown_per_combo = []
rbp_names_list = []

for idx, row in top_combo_df.iterrows():
    combo = row['v_combo']
    rbp_sub = df_merged[df_merged['v_combo'] == combo]
    rbp_names_list.append(", ".join(rbp_sub['RBP'].tolist()))

    fc_counts = rbp_sub[func_classes].isin(['+', '✓', 1, True]).sum()
    all_func_breakdown_per_combo.append(fc_counts.to_dict())

df_all_func_breakdown = pd.DataFrame(all_func_breakdown_per_combo)



# Plot layout and customisation logic
sns.set_style("ticks")
fig, axs = plt.subplots(
    nrows=2, ncols=2,
    figsize=(13.0, 7.0),
    gridspec_kw={
        'width_ratios': [2.2, 8.0],
        'height_ratios': [2.2, 1.8],
        'hspace': 0.05,
        'wspace': 0.06
    }
)

ax_dummy = axs[0, 0]
ax_bar = axs[0, 1]
ax_set_size = axs[1, 0]
ax_grid = axs[1, 1]

ax_dummy.axis('off')

x_limit_range = (-0.6, num_intersections - 0.4)

# Top-right panel: Stacked vertical bar chart for functions
bottoms = np.zeros(num_intersections)
for fc in func_classes:
    values = df_all_func_breakdown[fc].values
    if values.sum() > 0:
        ax_bar.bar(
            x_coordinates, values, bottom=bottoms, color=func_palette[fc],
            width=0.70, label=fc, zorder=3, edgecolor='white', linewidth=0.5
        )
        bottoms += values

sns.despine(ax=ax_bar, top=True, right=True, left=False, bottom=True)
ax_bar.grid(axis='y', linestyle='--', linewidth=0.5, color='#e5e8e8', zorder=0)
ax_bar.tick_params(axis='y', labelsize=12)

ax_bar.set_ylabel("Functional classes for RBPs", fontsize=12, fontweight='normal')
ax_bar.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax_bar.set_xlim(x_limit_range)

max_bar_height = max(bottoms) if len(bottoms) > 0 else 1
ax_bar.set_ylim(0, max_bar_height + 2.0)
ax_bar.yaxis.set_visible(False)
ax_bar.tick_params(axis='y', left=False, labelleft=False)
ax_bar.spines['left'].set_visible(False)
ax_bar.spines['right'].set_visible(False)
ax_bar.spines['top'].set_visible(False)

# Annotate RBP gene symbols above stacked bars
for x, rbp_str, total_h in zip(x_coordinates, rbp_names_list, bottoms):
    ax_bar.text(
        x, total_h + 0.15, rbp_str,
        ha='center', va='bottom', fontsize=12,
        rotation=90, color='black', fontweight='normal', zorder=5
    )

ax_bar.legend(
    title="RBP function", bbox_to_anchor=(1.02, 1), loc='upper left',
    frameon=False, fontsize=10, title_fontsize=10
)

# Bottom-left panel: Horizontal bar chart for virus set sizes
y_coordinates = np.arange(len(active_v_sets))
ax_set_size.barh(
    y_coordinates, total_v_sizes, color=left_set_bar_color,
    height=0.68, alpha=0.9, zorder=3
)
ax_set_size.invert_xaxis()
sns.despine(ax=ax_set_size, top=True, right=True, left=False, bottom=False)

ax_set_size.tick_params(axis='x', labelsize=12)
ax_set_size.set_xlabel("Total RBPs per virus family", fontsize=12, fontweight='normal')

ax_set_size.set_ylim(-0.5, len(active_v_sets) - 0.5)
ax_set_size.set_yticks(y_coordinates)
ax_set_size.set_yticklabels(active_v_sets, fontsize=10, fontweight='normal', fontstyle='italic')

max_val_set = max(total_v_sizes) if total_v_sizes else 1
ax_set_size.set_xlim(left=max_val_set * 1.22, right=0)

# Bottom-right panel: Intersection dot matrix grid
for y_idx in range(len(active_v_sets)):
    bg_color = '#f8f9f9' if y_idx % 2 == 0 else '#ffffff'
    ax_grid.axhspan(y_idx - 0.5, y_idx + 0.5, facecolor=bg_color, zorder=0)

for x_idx, row in top_combo_df.iterrows():
    combo = row['v_combo']
    active_indices = [y_idx for y_idx, act in enumerate(combo) if act]

    if len(active_indices) > 1:
        ax_grid.plot(
            [x_idx, x_idx], [min(active_indices), max(active_indices)],
            color="#2471A3", linewidth=1.2, zorder=1
        )

    inactive_indices = [y_idx for y_idx, act in enumerate(combo) if not act]
    ax_grid.scatter([x_idx] * len(inactive_indices), inactive_indices, color='#e5e8e8', s=25, zorder=2)
    ax_grid.scatter([x_idx] * len(active_indices), active_indices, color="#2471A3", s=25, zorder=3, edgecolors='none')

sns.despine(ax=ax_grid, top=True, right=True, left=True, bottom=False)

ax_grid.set_ylim(-0.5, len(active_v_sets) - 0.5)
ax_grid.set_yticks(y_coordinates)
ax_grid.set_yticklabels([])
ax_grid.tick_params(axis='y', length=0)

ax_grid.set_xlim(x_limit_range)
ax_grid.set_xticks(x_coordinates)
ax_grid.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

ax_grid.set_xlabel("\nShared RBP intersections", fontsize=12, fontweight='normal')



plt.show()
plt.close()
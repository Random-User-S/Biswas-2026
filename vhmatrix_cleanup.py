# This script does the following:
# 1. Mapping viruses into virus families (only animal viruses)
# 2. Combine into a matrix and outputs rbp_family_matrix.csv which is used for plotting heatmap


import pandas as pd
import re

# File paths
input_associations_path = "rbp_virus_associations.csv"
virus_patterns_path = "raw_to_clean_patterns.csv"
family_mapping_path = "clean_to_family.csv"


# Mapping logic:
#1. Load and pre-compile raw text to clean virus name regex
#2. Load clean virus name to virus family dictionary


def load_raw_to_clean_patterns(patterns_path=virus_patterns_path):
    df_patterns = pd.read_csv(patterns_path)
    return [
        (row["clean_name"], re.compile(row["regex_pattern"], re.IGNORECASE))
        for _, row in df_patterns.iterrows()
    ]


def load_clean_to_family_mapping(family_path=family_mapping_path):
    df_family = pd.read_csv(family_path)
    return dict(zip(df_family["clean_name"], df_family["virus_family"]))



compiled_patterns = load_raw_to_clean_patterns(virus_patterns_path)
family_lookup = load_clean_to_family_mapping(family_mapping_path)


def map_raw_virus_name(raw_name):
    if not isinstance(raw_name, str):
        return raw_name

    for clean_name, pattern in compiled_patterns:
        if pattern.search(raw_name):
            return clean_name

    return raw_name



# Load dataset
df = pd.read_csv(input_associations_path)

df["Virus_Clean"] = df["Virus"].apply(map_raw_virus_name)

df["Virus_Family"] = df["Virus_Clean"].map(family_lookup).fillna("Unclassified")



# Labeling and metadata
df["is_mapped"] = df["Virus_Family"] != "Unclassified"
df["Virus_Label"] = df.apply(
    lambda r: f'{r["Virus_Clean"]} ({r["Virus_Family"]})'
    if r["is_mapped"]
    else r["Virus"],
    axis=1,
)

# Separate into clean mapped entries and noisy unmapped entries
df_mapped = df[df["is_mapped"]].copy()
df_noisy = df[~df["is_mapped"]].copy()

# Export family matrix for downstream heatmap plotting
family_matrix = df_mapped.pivot_table(
    index="RBP",
    columns="Virus_Family",
    values="PMID",
    aggfunc="nunique",
    fill_value=0,
)
family_matrix.to_csv("rbp_family_matrix.csv")



# Export unmapped/noisy terms for manual pattern refinement
noisy_summary = (
    df_noisy.groupby("Virus")
    .agg(
        Total_Papers=("PMID", "nunique"),
        RBP_Count=("RBP", "nunique"),
        RBPs=("RBP", lambda x: ", ".join(sorted(x.unique()))),
    )
    .sort_values("Total_Papers", ascending=False)
    .reset_index()
)
noisy_summary.to_csv("rbp_virus_noisy.csv", index=False)


print(
    f"Mapped {len(df_mapped)} records across {df_mapped['Virus_Family'].nunique()} families."
)
print(f"Saved family matrix to 'rbp_family_matrix.csv'.")
print(f"Saved {len(noisy_summary)} unmapped terms to 'rbp_virus_noisy.csv'.")
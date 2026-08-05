# This script makes the combined summary 'rbp_master_literature_matrix.xlsx'

import re
import numpy as np
import pandas as pd



# File paths
assoc_csv_path = "rbp_virus_associations.csv"
names_csv_path = "rbp_full_names.csv"
patterns_csv_path = "raw_to_clean_patterns.csv"
family_csv_path = "clean_to_family.csv"
output_excel_path = "rbp_master_literature_matrix.xlsx"



def load_raw_to_clean_patterns(patterns_path=patterns_csv_path):
    df_patterns = pd.read_csv(patterns_path)
    return [
        (row["clean_name"], re.compile(row["regex_pattern"], re.IGNORECASE))
        for _, row in df_patterns.iterrows()
    ]


def load_clean_to_family_mapping(family_path=family_csv_path):
    df_family = pd.read_csv(family_path)
    family_dict = dict(zip(df_family["clean_name"], df_family["virus_family"]))
    all_families = sorted(df_family["virus_family"].unique().tolist())
    return family_dict, all_families


def map_raw_virus_name(raw_name, compiled_patterns):
    if not isinstance(raw_name, str):
        return raw_name

    # Pattern order in CSV
    for clean_name, pattern in compiled_patterns:
        if pattern.search(raw_name):
            return clean_name

    return raw_name




compiled_patterns = load_raw_to_clean_patterns(patterns_csv_path)
family_lookup, family_columns = load_clean_to_family_mapping(family_csv_path)

df_assoc = pd.read_csv(assoc_csv_path)
df_names = pd.read_csv(names_csv_path)

df_assoc["Virus_Clean"] = df_assoc["Virus"].apply(
    lambda x: map_raw_virus_name(x, compiled_patterns)
)
df_assoc["Virus_Family"] = df_assoc["Virus_Clean"].map(family_lookup).fillna("Unclassified")

# Safety check: Exclude unmapped terms before computing family-level metrics
df_human = df_assoc[df_assoc["Virus_Family"] != "Unclassified"].copy()

raw_matrix = df_human.pivot_table(
    index="RBP",
    columns="Virus_Family",
    values="PMID",
    aggfunc="nunique",
    fill_value=0,
)

# log1p transformation log(1 + x) prevents log(0) errors and compresses scale for a nice plot
log_matrix_families = np.log1p(raw_matrix).round(4)

rbp_summary = (
    df_human.groupby("RBP")
    .agg(
        Total_Papers_All_Families=("PMID", "nunique"),
        Specific_Viruses_Found=("Virus_Clean", lambda v: ", ".join(sorted(v.unique()))),
        PMID_List=("PMID", lambda p: ", ".join(map(str, sorted(p.unique())))),
    )
    .reset_index()
)

# Left-join against master gene frame ensures candidate RBPs with 0 hits remain in the final matrix
master_combined = pd.merge(
    df_names.drop(columns=["Rank"], errors="ignore"),
    rbp_summary,
    left_on="Query_RBP",
    right_on="RBP",
    how="left",
).drop(columns=["RBP"])

master_combined["Total_Papers_All_Families"] = (
    master_combined["Total_Papers_All_Families"].fillna(0).astype(int)
)
master_combined["Specific_Viruses_Found"] = master_combined["Specific_Viruses_Found"].fillna("None")
master_combined["PMID_List"] = master_combined["PMID_List"].fillna("N/A")

master_combined["Log1p_Total_Papers"] = np.log1p(
    master_combined["Total_Papers_All_Families"]
).round(4)

master_combined = pd.merge(
    master_combined,
    log_matrix_families,
    left_on="Query_RBP",
    right_index=True,
    how="left",
)


# Enforces matrix completeness for 0 hits
for col in family_columns:
    if col not in master_combined.columns:
        master_combined[col] = 0.0

master_combined[family_columns] = master_combined[family_columns].fillna(0.0)

master_combined = master_combined.sort_values(
    by="Total_Papers_All_Families", ascending=False
).reset_index(drop=True)

master_combined.insert(0, "Master_Rank", master_combined.index + 1)

master_combined.to_excel(output_excel_path, index=False)
print(f"\nSaved master literature matrix -> '{output_excel_path}'")
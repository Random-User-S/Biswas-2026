# This script maps the different processes into common functional classes
# The output is used for plotting the upset plot using the script 'upset_plot.py'



import re
import pandas as pd



# File paths
input_csv_path = "rbp_biological_processes.csv"
patterns_csv_path = "functional_patterns.csv"
output_csv_path = "rbp_functional_classes_matrix.csv" # This has the compiled regex patterns from the original csv


# Total evidence codes in analysed csv:
"""
IDA             
IEA             
ISS             
IMP             
IBA             
TAS             
NAS             
IC              
HDA             
IGI             
"""

# However, the other codes were excluded
allowed_evidence_codes = [
    "IDA",
    "IGI",
    "HDA",
    "ISS",
    "IBA",
    "IC",
    "TAS",
]

def load_functional_patterns(patterns_path=patterns_csv_path):

    df_patterns = pd.read_csv(patterns_path)
    return dict(zip(df_patterns["functional_class"], df_patterns["regex_pattern"]))



def classify_rbp_functions(
    input_csv=input_csv_path,
    patterns_csv=patterns_csv_path,
    output_csv=output_csv_path,
):
    # Load the regex functional patterns
    functional_patterns = load_functional_patterns(patterns_csv)

    df_go = pd.read_csv(input_csv)
    all_rbps = df_go["RBP_Symbol"].dropna().unique()

    # Filter out IEA and AES
    df_filtered = df_go[
        df_go["Evidence_Code"]
        .astype(str)
        .str.strip()
        .str.upper()
        .isin(allowed_evidence_codes)
    ].copy()

    # Aggregate all GO terms into a single lowercase text block per RBP
    rbp_grouped = (
        df_filtered.groupby("RBP_Symbol")["GO_Name"]
        .apply(lambda terms: " ".join(str(t) for t in terms).lower())
        .to_dict()
    )

    matrix_rows = []

    # Map GO terms to high-level functional classes
    for symbol in all_rbps:
        all_go_text = rbp_grouped.get(symbol, "")
        rbp_entry = {"RBP_Symbol": symbol}

        for f_class, pattern in functional_patterns.items():
            if all_go_text and re.search(pattern, all_go_text, re.IGNORECASE):
                rbp_entry[f_class] = "+"
            else:
                rbp_entry[f_class] = ""

        matrix_rows.append(rbp_entry)

    df_matrix = pd.DataFrame(matrix_rows)
    df_matrix.to_csv(output_csv, index=False)

    print(f"Output saved as '{output_csv}'")
    return df_matrix




df_matrix = classify_rbp_functions()
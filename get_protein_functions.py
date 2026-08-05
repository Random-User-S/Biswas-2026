# This script gets the RBP functions of all RBP Uniprot IDs
# The output is used for mapping the processes into common categories using the script 'GO_process_mapping.py'
# We query [https://rest.uniprot.org/uniprotkb/](https://rest.uniprot.org/uniprotkb/){accession}.json
# and get the json objects; then we extract the GO terms from the json objects


import json
import pandas as pd
import requests



# File paths
mapping_df = pd.read_csv("rbp_uniprot_mapping.csv")
output_csv = "rbp_biological_processes.csv"

# Read the input csv
rbp_symbols = (
    mapping_df["UniProtKB_CURIE"]
    .dropna()
    .drop_duplicates()
    .tolist()
)
print(f"Found {len(rbp_symbols)} UniProt IDs.")
print(rbp_symbols[:5])


# Get all functions using the UniProt REST API > json format
def get_rbp_bp_uniprot(gene_list):
    results = {}

    with requests.Session() as session:
        for curie_id in gene_list:
            accession = curie_id.replace("UniProtKB:", "").strip()
            url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"

            response = session.get(url, timeout=10)

            if response.status_code == 200:
                data = response.json()
                annotations = []
                seen_go_ids = set()

                for xref in data.get("uniProtKBCrossReferences", []):
                    if xref.get("database") == "GO":
                        go_id = xref.get("id")
                        properties = {
                            p.get("key"): p.get("value")
                            for p in xref.get("properties", [])
                        }
                        term_str = properties.get("GoTerm", "")

                        # Biological Process terms start with 'P:'
                        if term_str.startswith("P:") and go_id not in seen_go_ids:
                            seen_go_ids.add(go_id)
                            annotations.append({
                                "go_id": go_id,
                                "go_name": term_str[2:],  # Strip "P:"
                                "aspect": "Biological Process",
                                "evidence_code": properties.get(
                                    "GoEvidenceType", ""
                                ).split(":")[0],
                            })

                results[curie_id] = annotations
            else:
                results[curie_id] = []

    return results


# Json > table > csv
def export_go_results_to_csv(
    protein_functions,
    mapping_df,
    output_csv="rbp_biological_processes.csv",
):
    curie_to_symbol = dict(
        zip(mapping_df["UniProtKB_CURIE"], mapping_df["RBP_Symbol"])
    )

    rows = []

    for curie_id, annotations in protein_functions.items():
        gene_symbol = curie_to_symbol.get(curie_id, "N/A")

        if not annotations:
            rows.append({
                "RBP_Symbol": gene_symbol,
                "Mapped_UniProt_ID": curie_id,
                "GO_ID": "N/A",
                "GO_Name": "No Biological Process terms found",
                "Aspect": "Biological Process",
                "Evidence_Code": "N/A",
            })
        else:
            for ann in annotations:
                rows.append({
                    "RBP_Symbol": gene_symbol,
                    "Mapped_UniProt_ID": curie_id,
                    "GO_ID": ann.get("go_id"),
                    "GO_Name": ann.get("go_name"),
                    "Aspect": ann.get("aspect"),
                    "Evidence_Code": ann.get("evidence_code"),
                })

    df_results = pd.DataFrame(rows)
    df_results.to_csv(output_csv, index=False)

    print(
        f"Saved {len(df_results)} GO annotation records to "
        f"'{output_csv}'"
    )

    return df_results




# Main
print("Fetching GO Biological Processes from UniProt...")
protein_functions = get_rbp_bp_uniprot(rbp_symbols)

print("\nJSON Summary Sample:")
sample_keys = list(protein_functions.keys())[:2]
sample_data = {k: protein_functions[k] for k in sample_keys}
print(json.dumps(sample_data, indent=2))

print("\nOutput saved to csv:")
results_df = export_go_results_to_csv(
    protein_functions,
    mapping_df,
    output_csv=output_csv,
)
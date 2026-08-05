# This script convert all gene names to Uniprot IDs in the format: "UniProtKB:Q9NQ94"
# This output will be used to query and fetch the functions of RBPs from Uniprot using 'get_protein_functions.py'


import pandas as pd
import requests

# Input RBP gene list
rbp_df = pd.read_csv("rbp_list.csv")
my_rbp_symbols = (
    rbp_df["RBP_Symbol"]
    .dropna()
    .drop_duplicates()
    .tolist()
)


# Converts the gene names to Uniprot accessions
# Queries MyGene.info API to map gene symbols to UniProt accessions

def map_rbps_to_uniprot(gene_list, species="human"):

  url = "https://mygene.info/v3/query"

  # These are not working automatically
  MANUAL_LOOKUP = {
      "FMR1": "Q06787",
      "EEF1A": "P68104",  # EEF1A1
      "PTB3": "P26599",  # PTBP1
  }

  payload = {
      "q": ",".join(gene_list),
      "scopes": "symbol,alias,ensembl.gene,entrezgene",
      "fields": "uniprot",
      "species": species,
      "size": 1000,
  }

  response = requests.post(url, data=payload)
  mapping_results = []
  mapped_symbols = set()

  if response.status_code == 200:
    for item in response.json():
      query = item.get("query")
      if not query or query in mapped_symbols:
        continue

      uniprot_data = item.get("uniprot", {})
      primary_id = None

      # Check nested dictionary structure
      if isinstance(uniprot_data, dict):
        sp = uniprot_data.get("Swiss-Prot")
        trembl = uniprot_data.get("TrEMBL")
        entry = sp or trembl
        if isinstance(entry, list) and len(entry) > 0:
          primary_id = entry[0]
        elif isinstance(entry, str):
          primary_id = entry

      # Check flat string structure
      elif isinstance(uniprot_data, str):
        primary_id = uniprot_data

      if not primary_id and query in MANUAL_LOOKUP:
        primary_id = MANUAL_LOOKUP[query]

      if primary_id:
        mapping_results.append({
            "RBP_Symbol": query,
            "UniProt_Accession": primary_id,
            "UniProtKB_CURIE": f"UniProtKB:{primary_id}",
        })
        mapped_symbols.add(query)

  # Check unmapped queries against manual lookup
  for symbol in gene_list:
    if symbol not in mapped_symbols:
      pid = MANUAL_LOOKUP.get(symbol)
      if pid:
        mapping_results.append({
            "RBP_Symbol": symbol,
            "UniProt_Accession": pid,
            "UniProtKB_CURIE": f"UniProtKB:{pid}",
        })
        mapped_symbols.add(symbol)
      else:
        print(f"Failed: '{symbol}'")

  return pd.DataFrame(mapping_results)



# Main
df_mapped = map_rbps_to_uniprot(my_rbp_symbols)



csv_filename = "rbp_uniprot_mapping.csv"
df_mapped.to_csv(csv_filename, index=False)
print(f"Output saved to '{csv_filename}'\n")


# Format into quoted, comma-separated string
quoted_curies = [f'"{curie}"' for curie in df_mapped["UniProtKB_CURIE"]]
formatted_string = ", ".join(quoted_curies)

print("Mapped UniProt IDs:")
print(formatted_string)
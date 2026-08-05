# This script does the following:
# 1. Input RBP list
# 2. Get research papers from Pubmed using custom a query
# 3. Retrieve PubTator3 annotations for the papers to identify annotate the co-occurences of RBPs and viruses
# 4. Save data matrices as csv files


import pandas as pd
import requests
import time


max_pmids = 10000
batch_size = 100


base_esearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
base_efetch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
base_pubtator = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"


rbp_df = pd.read_csv("rbp_list.csv")
rbps = rbp_df["RBP_Symbol"].dropna().tolist()

# For caching logic
taxid_cache = {}


def is_virus(name, tax_id):
    name_lower = name.lower()

    if name_lower in non_virus_terms:
        return False

    tokens = set(name_lower.replace("-", " ").split())
    if any(kw in name_lower for kw in virus_terms) or any(
            t in virus_terms for t in tokens
    ):
        return True

    # If tax_id is ok and in cache, then get it from cache. Otherwise query NCBI.
    if tax_id and tax_id.isdigit():
        if tax_id not in taxid_cache:
            base_esummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            params = {"db": "taxonomy", "id": tax_id, "retmode": "json"}
            r = requests.get(base_esummary, params=params, timeout=10)

            res = r.json().get("result")
            lineage = res.get(tax_id, {}).get("lineage", "") if isinstance(res, dict) else ""

            taxid_cache[tax_id] = "Viruses" in lineage
            time.sleep(0.5)
        return taxid_cache[tax_id]

    return False

# Queryind PubMed
def get_pmids(rbp):
    query = f'("{rbp}"[Title/Abstract]) AND ("Viruses"[MeSH Terms] OR "virus"[Title/Abstract] OR "viral"[Title/Abstract]) NOT "Review"[Publication Type]'
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": max_pmids,
    }
    r = requests.get(base_esearch, params=params, timeout=15)
    r.raise_for_status()
    return r.json()["esearchresult"]["idlist"]



def load_terms(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        return {line.strip().lower() for line in f if line.strip()}
        
# Use them for PubTator3
non_virus_terms = load_terms("non_virus_terms.txt")
virus_terms = load_terms("virus_terms.txt")

def get_virus_mentions(pmids):
    if not pmids:
        return []
    all_viruses = []

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        r = requests.get(
            base_pubtator,
            params={"pmids": ",".join(batch), "concepts": "Species"},
            timeout=30,
        )
        if r.status_code != 200:
            continue
        """
        Json structure:
        doc (Publication)
          └── passages (e.g., Passage 0 = Title, Passage 1 = Abstract)
                └── annotations (Tagged entities: Genes, Species, etc.)
        Collectively append to the list as "pmid": pmid, "name": name, "tax_id": tax_id
        """
        data = r.json()
        docs = data if isinstance(data, list) else data.get("PubTator3", [])
        for doc in docs:
            pmid = doc.get("id", "?")
            for passage in doc.get("passages", []):
                for ann in passage.get("annotations", []):
                    info = ann.get("infons", {})
                    if info.get("type") != "Species":
                        continue
                    name = ann.get("text", "").strip()
                    tax_id = info.get("identifier", "")
                    if name and is_virus(name, tax_id):
                        all_viruses.append(
                            {"pmid": pmid, "name": name, "tax_id": tax_id}
                        )

        time.sleep(0.5)

    return all_viruses




# Main
results = []

for rbp in rbps:
    pmids = get_pmids(rbp)
    viruses = get_virus_mentions(pmids)

    for v in viruses:
        results.append({"RBP": rbp, "PMID": v["pmid"], "Virus": v["name"]})

    time.sleep(0.5)




# Output rbp-virus associations
df = (
    pd.DataFrame(results)
    .drop_duplicates(subset=["RBP", "PMID", "Virus"])
    .sort_values(["RBP", "Virus"])
    .reset_index(drop=True)
)
print("Output rbp-virus associations saved to rbp_virus_associations.csv")
df.to_csv("rbp_virus_associations.csv", index=False)

# Output virus summary
summary = (
    df.groupby("Virus")
    .agg(
        Total_Papers=("PMID", "nunique"),
        RBP_Count=("RBP", "nunique"),
        RBPs=("RBP", lambda x: ", ".join(sorted(x.unique()))),
    )
    .sort_values("Total_Papers", ascending=False)
    .reset_index()
)
print("Output virus summary saved to virus_summary.csv")
summary.to_csv("virus_summary.csv", index=False)

# Output RBP-virus association matrix
matrix = df.pivot_table(
    index="RBP",
    columns="Virus",
    values="PMID",
    aggfunc="nunique",
    fill_value=0,
)
print("Output RBP-virus association matrix saved to rbp_virus_matrix.csv")
matrix.to_csv("rbp_virus_matrix.csv")
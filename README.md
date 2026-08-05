## Code used in "Host RNA-binding proteins as broad-spectrum targets for antiviral therapy"

This repository contains the Python scripts used to generate **Figure 1b (heatmap)** and **Figure 1c (UpSet plot)** for the manuscript.

---

### What is this?

The scripts get papers from PubMed and visualise the bibliographic associations of viruses and host RNA-binding proteins (RBPs). The annotated functions of the RBPs are also fetched from UniProt.

---

### Requirements

Before running the scripts, make sure you have:

- A computer
- Python **3.9** or later
- A copy of this repository
- The required Python packages in `requirements.txt`

---

### Workflow

You can run `controller.py`, which executes all these scripts sequentially. The run takes several minutes (~15 min) to complete.

Alternatively, you can run the scripts individually.

#### Figure 1b – Heatmap

Run the following scripts **in order**:

1. `get_viruses_rbps.py`
2. `vhmatrix_cleanup.py`
3. `heatmap_plot.py`



#### Figure 1c – UpSet plot

Run the following scripts **in order**:

1. `summary_combined.py`
2. `convert_protein_IDs.py`
3. `get_protein_functions.py`
4. `GO_process_mapping.py`
5. `upset_plot.py`



---

#### Notes:

1. Data was retrieved from PubMed and PubTator3 on July 13, 2026. The exact dataset generated on that date is supplied as `rbp_master_literature_matrix_v1.xlsx`.
2. Although this repository is configured for RBPs, `get_viruses_rbps.py` is a general-purpose literature mining script that can be used to retrieve associations between viruses and any host protein.
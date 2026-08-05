import subprocess
import sys



# List the scripts in the exact order you want them to run
scripts = ["get_viruses_rbps.py","vhmatrix_cleanup.py","heatmap_plot.py","summary_combined.py",
           "convert_protein_IDs.py","get_protein_functions.py","GO_process_mapping.py","upset_plot.py"
]

for script in scripts:
    print(f"--- Starting {script} ---")

    # sys.executable points to the current Python interpreter
    result = subprocess.run([sys.executable, script], capture_output=False, text=True)

    # Check if the script executed successfully
    if result.returncode == 0:
        print(f"--- Finished {script} successfully ---\n")
    else:
        print(f"❌ Error: {script} failed with exit code {result.returncode}. Stopping master execution.")
        break  # Optional: stop the master script if one step fails
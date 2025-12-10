import pandas as pd
import os

# -------------------------------
# File paths
# -------------------------------
genes_file = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/ClinGen/Natural History/included_inheritance_genesymbols.txt"

mgi_file = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Raw Datasets/MGI/mouse_model_disease_clean_ds.txt"

output_folder = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/ClinGen/Natural History/Models"
os.makedirs(output_folder, exist_ok=True)

mgi_csv_output = os.path.join(output_folder, "mgi_models.csv")
mgi_txt_output = os.path.join(output_folder, "mgi_genesymbols.txt")
modeled_genes_txt = os.path.join(output_folder, "mgi_modeled_genesymbols.txt")

# -------------------------------
# Check that files exist
# -------------------------------
for file_path in [genes_file, mgi_file]:
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

# -------------------------------
# Load included genes
# -------------------------------
with open(genes_file, "r") as f:
    included_genes = {line.strip() for line in f if line.strip() != ""}

# -------------------------------
# Load MGI dataset
# -------------------------------
mgi_df = pd.read_csv(mgi_file, sep="\t", encoding="utf-8")

# Strip leading/trailing spaces from HumanHomologs
mgi_df["HumanHomologs"] = mgi_df["HumanHomologs"].astype(str).str.strip()

# -------------------------------
# Cross-reference genes
# -------------------------------
found_df = mgi_df[mgi_df["HumanHomologs"].isin(included_genes)].copy()

# Ensure MouseModels is numeric
found_df["MouseModels"] = pd.to_numeric(found_df["MouseModels"], errors="coerce").fillna(0)

# Genes with at least 1 mouse model
modeled_df = found_df[found_df["MouseModels"] >= 1]

# -------------------------------
# Save CSV with selected columns
# -------------------------------
selected_columns = ["DiseaseTerm", "MouseHomologs", "HumanHomologs", "MouseModels", "HomologySource"]
modeled_df.to_csv(mgi_csv_output, columns=selected_columns, index=False)

# -------------------------------
# Save TXT with all found genes
# -------------------------------
found_genes = sorted(found_df["HumanHomologs"].unique())
with open(mgi_txt_output, "w") as f:
    for gene in found_genes:
        f.write(gene + "\n")

# -------------------------------
# Save TXT with genes with at least 1 mouse model
# -------------------------------
modeled_genes = sorted(modeled_df["HumanHomologs"].unique())
with open(modeled_genes_txt, "w") as f:
    for gene in modeled_genes:
        f.write(gene + "\n")

# -------------------------------
# Terminal output (visual summary)
# -------------------------------
total_found_genes = found_df["HumanHomologs"].nunique()
genes_with_models = modeled_df["HumanHomologs"].nunique()
total_rows = len(modeled_df)

print(f"✅ Found {total_found_genes} unique gene symbols in MGI dataset after flexible parsing.")
print(f"💾 Saved filtered dataset to: {mgi_csv_output}")
print(f"💾 Saved gene symbols to: {mgi_txt_output}")
print(f"💾 Saved modeled gene symbols to: {modeled_genes_txt}\n")

print("=== SUMMARY ===")
print(f"Filter gene symbols provided: {len(included_genes)}")
print(f"Found in MGI dataset: {total_found_genes}")
print(f"Genes with at least 1 mouse model: {genes_with_models}")
print(f"Final rows in output CSV: {total_rows}")
print("✅ Process complete.")

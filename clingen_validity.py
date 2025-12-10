import pandas as pd
import os

# =====================
# FILE PATHS
# =====================
CLINGEN_FILE = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/ClinGen/clingen_gene_disease_validity_ds.csv"
ORPHANET_FILE = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/kept_disorders.csv"
OUTPUT_FOLDER = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/ClinGen"

os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# =====================
# LOAD CLINGEN CSV
# =====================
# Step 1: Load without header to find the actual header row
raw_df = pd.read_csv(CLINGEN_FILE, header=None)
header_row_index = None
for i, row in raw_df.iterrows():
    if "GENE SYMBOL" in row.values:
        header_row_index = i
        break

if header_row_index is None:
    raise ValueError("Could not find header row containing 'GENE SYMBOL'")

# Step 2: Reload CSV using the correct header
clingen_df = pd.read_csv(CLINGEN_FILE, header=header_row_index)
clingen_df.columns = clingen_df.columns.str.strip()  # remove extra whitespace
clingen_df["GENE SYMBOL"] = clingen_df["GENE SYMBOL"].str.strip()

# =====================
# LOAD ORPHANET CSV
# =====================
orphanet_df = pd.read_csv(ORPHANET_FILE)
orphanet_df["GeneSymbol"] = orphanet_df["GeneSymbol"].astype(str).str.strip()

print(f"Loaded {len(orphanet_df)} Orphanet disorders")
unique_orphanet_genes = set(orphanet_df["GeneSymbol"])
print(f"Unique gene symbols from Orphanet: {len(unique_orphanet_genes)}")

# =====================
# MATCH ORPHANET GENES TO CLINGEN
# =====================
matched_df = clingen_df[clingen_df["GENE SYMBOL"].isin(unique_orphanet_genes)].copy()
matched_genes = set(matched_df["GENE SYMBOL"])
unmatched_genes = unique_orphanet_genes - matched_genes

print(f"Orphanet disorders with at least one ClinGen match: {len(matched_df)}")

# Save unmatched genes
no_match_genes_file = os.path.join(OUTPUT_FOLDER, "no_match_genes.txt")
with open(no_match_genes_file, "w") as f:
    for gene in sorted(unmatched_genes):
        f.write(f"{gene}\n")

# Save unmatched orphanet disorders
orphanet_no_match_df = orphanet_df[orphanet_df["GeneSymbol"].isin(unmatched_genes)]
orphanet_no_match_csv = os.path.join(OUTPUT_FOLDER, "orphanet_no_match.csv")
orphanet_no_match_df.to_csv(orphanet_no_match_csv, index=False)

# =====================
# DEBUG: ClinGen classification counts
# =====================
print("ClinGen classification counts:")
print(matched_df["CLASSIFICATION"].value_counts())

# =====================
# FILTER DEFINITIVE
# =====================
definitive_df = matched_df[matched_df["CLASSIFICATION"].str.upper() == "DEFINITIVE"].copy()
definitive_genes = set(definitive_df["GENE SYMBOL"])

print(f"Number of unique gene symbols with Definitive status: {len(definitive_genes)}")

# Cross reference with Orphanet to get definitive OrphaCodes
definitive_orphanet_df = orphanet_df[orphanet_df["GeneSymbol"].isin(definitive_genes)]
definitive_orphacodes = set(definitive_orphanet_df["OrphaCode"])
print(f"Number of Orphanet entities with definitive ClinGen classification: {len(definitive_orphacodes)}")

# =====================
# SAVE OUTPUT FILES
# =====================
# All matches
matched_csv = os.path.join(OUTPUT_FOLDER, "clingen_all_matches.csv")
matched_df.to_csv(matched_csv, index=False)

# Definitive matches
definitive_csv = os.path.join(OUTPUT_FOLDER, "clingen_definitive.csv")
definitive_df.to_csv(definitive_csv, index=False)

# Definitive gene symbols
definitive_genes_txt = os.path.join(OUTPUT_FOLDER, "clingen_definitive_genesymbols.txt")
with open(definitive_genes_txt, "w") as f:
    for gene in sorted(definitive_genes):
        f.write(f"{gene}\n")

# Definitive OrphaCodes
definitive_orphacodes_txt = os.path.join(OUTPUT_FOLDER, "definitive_orphacodes.txt")
with open(definitive_orphacodes_txt, "w") as f:
    for oc in sorted(definitive_orphacodes):
        f.write(f"{oc}\n")

print(f"No match genes TXT saved: {no_match_genes_file}")
print(f"No match Orphanet CSV saved: {orphanet_no_match_csv}")
print(f"All ClinGen matches saved: {matched_csv}")
print(f"Definitive ClinGen matches saved: {definitive_csv}")
print(f"Definitive gene symbols TXT saved: {definitive_genes_txt}")
print(f"Definitive OrphaCodes TXT saved: {definitive_orphacodes_txt}")
print("✅ ClinGen validity filtering pipeline complete!")

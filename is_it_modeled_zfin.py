#!/usr/bin/env python3
import os
import pandas as pd

# === FILE PATHS ===
BASE_DIR = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders"

ZFIN_DIR = os.path.join(BASE_DIR, "Raw Datasets/ZFIN")
CLINGEN_NH_DIR = os.path.join(BASE_DIR, "Gene Association/Monogenic Association/Disorder Classification/ClinGen/Natural History")
OUTPUT_DIR = os.path.join(CLINGEN_NH_DIR, "Models")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Input files
GENE_DISEASE_FILE = os.path.join(ZFIN_DIR, "gene2DiseaseViaOrthology_ds.txt")
FISH_MODEL_FILE = os.path.join(ZFIN_DIR, "fish_model_disease_ds.txt")
GENE_FILTER_FILE = os.path.join(CLINGEN_NH_DIR, "included_inheritance_genesymbols.txt")

# Output files
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "zfin_models.csv")
OUTPUT_TXT = os.path.join(OUTPUT_DIR, "zfin_genesymbols.txt")

# === STEP 1: Load filter genes ===
with open(GENE_FILTER_FILE, "r") as f:
    filter_genes = {line.strip() for line in f if line.strip()}

print(f"Loaded {len(filter_genes)} gene symbols from filter file.")

# === STEP 2: Load ZFIN gene→disease orthology data ===
df_gene_disease = pd.read_csv(GENE_DISEASE_FILE, sep="\t", dtype=str)
df_gene_disease.columns = df_gene_disease.columns.str.strip()

print(f"Loaded {len(df_gene_disease)} rows from gene2DiseaseViaOrthology_ds.txt")

# === STEP 3: Filter for relevant human orthologs ===
filtered_genes = df_gene_disease[df_gene_disease["HumanOrthologSymbol"].isin(filter_genes)].copy()
found_genes = filtered_genes["HumanOrthologSymbol"].nunique()

print(f"✅ Found {found_genes} unique gene symbols in ZFIN gene→disease dataset.")
print(f"Total matching rows kept: {len(filtered_genes)}")

# === STEP 4: Extract DOTermIDs to cross-reference ===
doterm_ids = filtered_genes["DOTermID"].dropna().unique().tolist()

# === STEP 5: Load ZFIN fish model dataset ===
df_fish_model = pd.read_csv(FISH_MODEL_FILE, sep="\t", dtype=str)
df_fish_model.columns = df_fish_model.columns.str.strip()

print(f"Loaded {len(df_fish_model)} rows from fish_model_disease_ds.txt")

# === STEP 6: Filter fish model data by DOTermID ===
filtered_models = df_fish_model[df_fish_model["DOTermID"].isin(doterm_ids)].copy()
print(f"✅ Found {len(filtered_models)} matching rows in fish model dataset.")
print(f"Unique DOTermIDs with models: {filtered_models['DOTermID'].nunique()}")

# === STEP 7: Merge the two datasets on DOTermID ===
merged = pd.merge(
    filtered_genes,
    filtered_models,
    on="DOTermID",
    how="inner",
    suffixes=("_orthology", "_model")
)

print(f"✅ Total merged rows: {len(merged)}")

# === STEP 8: Select required columns ===
expected_columns = [
    "HumanOrthologSymbol",
    "DOTermID",
    "DOTermName_orthology",  # from first dataset
    "OMIMTermName",
    "OMIMID",
    "is_a_model",            # from fish model dataset
    "PubMedID"
]

available_columns = [col for col in expected_columns if col in merged.columns]
missing_columns = set(expected_columns) - set(available_columns)
if missing_columns:
    print(f"⚠️ Warning: missing columns: {', '.join(missing_columns)}")

merged = merged[available_columns]

# === STEP 9: Save output CSV ===
merged.to_csv(OUTPUT_CSV, index=False)
print(f"💾 Saved merged ZFIN model dataset to: {OUTPUT_CSV}")

# === STEP 10: Save unique gene symbols ===
unique_genes = sorted(merged["HumanOrthologSymbol"].dropna().unique())
with open(OUTPUT_TXT, "w", encoding="utf-8") as f:
    for g in unique_genes:
        f.write(g + "\n")

print(f"💾 Saved unique gene symbols to: {OUTPUT_TXT}")
print(f"✅ Total unique modeled genes: {len(unique_genes)}")

# === FINAL SUMMARY ===
print("\n=== SUMMARY ===")
print(f"Filter gene symbols provided: {len(filter_genes)}")
print(f"Found in ZFIN orthology dataset: {found_genes}")
print(f"Found in ZFIN model dataset (DOTermID match): {filtered_models['DOTermID'].nunique()}")
print(f"Final merged rows: {len(merged)}")
print(f"Unique gene symbols modeled: {len(unique_genes)}")
print("✅ Process complete.")

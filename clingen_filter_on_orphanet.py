#!/usr/bin/env python3

import pandas as pd

# === Input paths ===
orphacodes_txt = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/ClinGen/definitive_orphacodes.txt"
input_csv = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/kept_disorders.csv"

# === Output path ===
output_csv = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/ClinGen/definitive_disorders.csv"

# === Load OrphaCodes from TXT file ===
with open(orphacodes_txt, "r") as f:
    filter_orphacodes = {line.strip() for line in f if line.strip()}  # use set for speed

print(f"Loaded {len(filter_orphacodes)} OrphaCodes from filter list")

# === Load input CSV ===
df = pd.read_csv(input_csv, dtype={"OrphaCode": str})  # keep OrphaCode as string
print(f"Loaded {len(df)} rows from input CSV")

# === Filter rows ===
filtered_df = df[df["OrphaCode"].astype(str).isin(filter_orphacodes)]
print(f"Filtered dataset contains {len(filtered_df)} rows")

# === Save output CSV ===
filtered_df.to_csv(output_csv, index=False)
print(f"Filtered file saved to: {output_csv}")

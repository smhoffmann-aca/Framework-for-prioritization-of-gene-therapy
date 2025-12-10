import xml.etree.ElementTree as ET
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# =====================
# CONFIGURATION
# =====================
OUTPUT_DIR = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Merge Neuro, IME, Genetic"

DATASETS = {
    "neurological": "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Raw Datasets/neurological_disorders_ds.xml",
    "genetic": "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Raw Datasets/genetic_diseases_ds.xml",
    "metabolic": "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Raw Datasets/inborn_errors_metabolism_ds.xml"
}

os.makedirs(OUTPUT_DIR, exist_ok=True)

# =====================
# FUNCTIONS
# =====================
def parse_xml(file_path):
    """Extract disorders with OrphaCode, Name, ExpertLink from an Orphanet XML file"""
    tree = ET.parse(file_path)
    root = tree.getroot()

    # DEBUG: show first few tags
    print(f"Inspecting {file_path} ... root tag = {root.tag}")
    for elem in root.iter():
        print("Tag:", elem.tag)
        break  # just show the first one

    records = []
    for disorder in root.findall(".//Disorder"):
        orpha = disorder.findtext("OrphaCode")
        name = disorder.findtext("Name[@lang='en']") or disorder.findtext("Name")
        link = disorder.findtext("ExpertLink[@lang='en']") or disorder.findtext("ExpertLink")
        
        if orpha:
            records.append({
                "OrphaCode": orpha.strip(),
                "Name": (name or "").strip(),
                "ExpertLink": (link or "").strip(),
                "SourceFile": os.path.basename(file_path)
            })

    print(f"Parsed {len(records)} disorders from {os.path.basename(file_path)}")
    return pd.DataFrame(records)


def save_tier(orpha_set, tier_name, all_data):
    """Save tier to CSV and TXT"""
    tier_df = all_data[all_data["OrphaCode"].isin(orpha_set)]
    csv_path = os.path.join(OUTPUT_DIR, f"{tier_name}.csv")
    txt_path = os.path.join(OUTPUT_DIR, f"{tier_name}.txt")

    tier_df.to_csv(csv_path, index=False)
    with open(txt_path, "w") as f:
        for oc in orpha_set:
            f.write(f"{oc}\n")

    print(f"Saved {tier_name}: {len(orpha_set)} disorders")


# =====================
# MAIN SCRIPT
# =====================
print("Parsing XML datasets...")

df_neuro = parse_xml(DATASETS["neurological"])
df_genetic = parse_xml(DATASETS["genetic"])
df_metabolic = parse_xml(DATASETS["metabolic"])

print(f"Neurological: {len(df_neuro)} disorders")
print(f"Genetic: {len(df_genetic)} disorders")
print(f"Metabolic: {len(df_metabolic)} disorders")

# Merge all datasets into a single master file
combined = pd.concat([df_neuro, df_genetic, df_metabolic], ignore_index=True)

# Normalize column names
combined.rename(columns={"ORPHACode": "OrphaCode"}, inplace=True)
combined.columns = combined.columns.str.strip()

# Keep only unique OrphaCodes, merging names/links if duplicated
combined = combined.groupby("OrphaCode").agg({
    "Name": "first",
    "ExpertLink": "first",
    "SourceFile": lambda x: ";".join(set(x))
}).reset_index()

# Save master CSV
master_csv = os.path.join(OUTPUT_DIR, "merged_master.csv")
combined.to_csv(master_csv, index=False)
print(f"Master dataset saved: {master_csv}")

# =====================
# TIERED LISTS
# =====================
set_neuro = set(df_neuro["OrphaCode"])
set_genetic = set(df_genetic["OrphaCode"])
set_metabolic = set(df_metabolic["OrphaCode"])

# Tier definitions
tier1 = set_neuro & set_genetic & set_metabolic
tier2 = (set_neuro & set_genetic) | (set_neuro & set_metabolic) | (set_genetic & set_metabolic)
tier2 = tier2 - tier1
tier3 = set_neuro | set_genetic | set_metabolic
tier3 = tier3 - (tier1 | tier2)

# Save tiered outputs
save_tier(tier1, "tier1_all_three", combined)
save_tier(tier2, "tier2_two_overlap", combined)
save_tier(tier3, "tier3_one_only", combined)

# =====================
# VENN DIAGRAM
# =====================
plt.figure(figsize=(8,8))
venn3(
    [set_neuro, set_genetic, set_metabolic],
    set_labels=("Neurological", "Genetic", "Metabolic")
)
venn_path = os.path.join(OUTPUT_DIR, "venn_diagram.png")
plt.savefig(venn_path)
print(f"Venn diagram saved: {venn_path}")

print("✅ Done!")

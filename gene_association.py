import xml.etree.ElementTree as ET
import csv
import os

# File paths
tier1_txt = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Merge Neuro, IME, Genetic/tier1_all_three.txt"
gene_assoc_xml = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Raw Datasets/genes_associated_ds.xml"
output_folder = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association"

os.makedirs(output_folder, exist_ok=True)

# Load Tier1 OrphaCodes
with open(tier1_txt, "r") as f:
    tier1_orphacodes = {line.strip() for line in f.readlines()}
print(f"Loaded {len(tier1_orphacodes)} Tier1 OrphaCodes from TXT")

# Parse XML
tree = ET.parse(gene_assoc_xml)
root = tree.getroot()

# Group containers
group1 = []  # Strong candidates
group2 = []  # Supplementary
group3 = []  # Excluded
group4_not_yet_assessed = []  # Not yet assessed

# Track matched codes
tier1_orphacodes_matched = set()
tier1_orphacodes_not_assessed = set()

# Counters
total_assessed = 0
total_not_assessed = 0
disorders_found = 0

# Group type definitions
group1_types = {
    "Disease-causing germline mutation(s) in",
    "Disease-causing germline mutation(s) (loss of function) in",
    "Disease-causing germline mutation(s) (gain of function) in"
}
group2_types = {
    "Disease-causing somatic mutation(s) in",
    "Modifying germline mutation in"
}
group3_types = {
    "Part of a fusion gene in",
    "Major susceptibility factor in",
    "Candidate gene tested in",
    "Biomarker tested in"
}

# Iterate disorders
for disorder in root.find(".//DisorderList"):
    orpha_code = disorder.findtext("OrphaCode", "").strip()
    if orpha_code not in tier1_orphacodes:
        continue
    tier1_orphacodes_matched.add(orpha_code)
    disorders_found += 1

    disorder_name = disorder.findtext("Name", "").strip()
    disorder_type_name = disorder.find("DisorderType/Name").text.strip() if disorder.find("DisorderType/Name") is not None else ""
    disorder_group_name = disorder.find("DisorderGroup/Name").text.strip() if disorder.find("DisorderGroup/Name") is not None else ""
    expert_link = disorder.findtext("ExpertLink", "").strip()

    # Count number of associations for this disorder
    gene_assoc_list = disorder.find("DisorderGeneAssociationList")
    assoc_count = int(gene_assoc_list.attrib.get("count", "0")) if gene_assoc_list is not None else 0
    if gene_assoc_list is None:
        continue

    for gene_assoc in gene_assoc_list.findall("DisorderGeneAssociation"):
        status = gene_assoc.findtext("DisorderGeneAssociationStatus/Name", "").strip()

        # Separate handling for Assessed vs Not yet assessed
        if status == "Assessed":
            total_assessed += 1
        elif status == "Not yet assessed":
            total_not_assessed += 1
            row = [orpha_code, disorder_name, disorder_group_name, disorder_type_name,
                   "", "", "", "", "", status, "", "", "", expert_link, assoc_count]
            group4_not_yet_assessed.append(row)
            tier1_orphacodes_not_assessed.add(orpha_code)
            continue  # don’t put into Groups 1–3
        else:
            continue

        assoc_type = gene_assoc.findtext("DisorderGeneAssociationType/Name", "").strip()

        # Gene info
        gene = gene_assoc.find("Gene")
        gene_symbol, gene_name, gene_type, gene_locus = "", "", "", ""
        omim_ref, uniprot_ref = "", ""
        source_of_validation = gene_assoc.findtext("SourceOfValidation", "").strip()

        if gene is not None:
            gene_symbol = gene.findtext("Symbol", "").strip()
            gene_name = gene.findtext("Name", "").strip()
            gene_type = gene.findtext("GeneType/Name", "").strip()

            locus_list = gene.find("LocusList")
            if locus_list is not None and int(locus_list.attrib.get("count", "0")) > 0:
                locus = locus_list.find("Locus/GeneLocus")
                if locus is not None:
                    gene_locus = locus.text.strip()

            ext_refs = gene.find("ExternalReferenceList")
            if ext_refs is not None:
                for ref in ext_refs.findall("ExternalReference"):
                    source = ref.findtext("Source", "").strip()
                    value = ref.findtext("Reference", "").strip()
                    if source == "OMIM":
                        omim_ref = value
                    if source in ["UNIPROTKB", "SwissProt"]:
                        uniprot_ref = value

        row = [
            orpha_code, disorder_name, disorder_group_name, disorder_type_name,
            gene_symbol, gene_name, gene_type, gene_locus,
            assoc_type, status, source_of_validation,
            omim_ref, uniprot_ref, expert_link, assoc_count
        ]

        # Assign to groups
        if assoc_type in group1_types:
            group1.append(row)
        elif assoc_type in group2_types:
            group2.append(row)
        elif assoc_type in group3_types:
            group3.append(row)

# Debug prints
print(f"Disorders in XML matching Tier1 OrphaCodes: {disorders_found}")
print(f"Total gene associations with status 'Assessed': {total_assessed}")
print(f"Total gene associations with status 'Not yet assessed': {total_not_assessed}")

# Save CSVs
headers = [
    "OrphaCode", "DisorderName", "DisorderGroup", "DisorderType",
    "GeneSymbol", "GeneName", "GeneType", "GeneLocus",
    "DisorderGeneAssociationType", "DisorderGeneAssociationStatus",
    "SourceOfValidation", "OMIM", "UniProt", "ExpertLink", "AssocCount"
]

group_files = [
    ("Group1_Strong.csv", group1),
    ("Group2_Supplementary.csv", group2),
    ("Group3_Excluded.csv", group3),
    ("Group4_NotYetAssessed.csv", group4_not_yet_assessed)
]

for fname, data in group_files:
    path = os.path.join(output_folder, fname)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(data)
    print(f"Saved {fname}: {len(data)} rows")

# Save TXT of matched OrphaCodes
matched_orphacodes_file = os.path.join(output_folder, "tier1_matched_orphacodes.txt")
with open(matched_orphacodes_file, "w", encoding="utf-8") as f:
    for code in sorted(tier1_orphacodes_matched):
        f.write(f"{code}\n")
print(f"Saved TXT of matched OrphaCodes: {matched_orphacodes_file}")

# Save TXT of not yet assessed OrphaCodes
not_assessed_file = os.path.join(output_folder, "tier1_not_yet_assessed_orphacodes.txt")
with open(not_assessed_file, "w", encoding="utf-8") as f:
    for code in sorted(tier1_orphacodes_not_assessed):
        f.write(f"{code}\n")
print(f"Saved TXT of not yet assessed OrphaCodes: {not_assessed_file}")

# NEW: Save unmatched OrphaCodes
unmatched_orphacodes = tier1_orphacodes - tier1_orphacodes_matched
unmatched_file = os.path.join(output_folder, "tier1_unmatched_orphacodes.txt")
with open(unmatched_file, "w", encoding="utf-8") as f:
    for code in sorted(unmatched_orphacodes):
        f.write(f"{code}\n")
print(f"Saved TXT of unmatched OrphaCodes: {unmatched_file}")
print(f"Total unmatched OrphaCodes: {len(unmatched_orphacodes)}")

print("✅ Finished gene association pipeline")

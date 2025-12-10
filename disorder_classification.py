import xml.etree.ElementTree as ET
import csv
import os

# File paths
monogenic_txt = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/monogenic_orphacodes.txt"
gene_assoc_xml = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Raw Datasets/genes_associated_ds.xml"
output_folder = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification"

os.makedirs(output_folder, exist_ok=True)

# Load monogenic OrphaCodes
with open(monogenic_txt, "r") as f:
    monogenic_orphacodes = {line.strip() for line in f.readlines()}
print(f"Loaded {len(monogenic_orphacodes)} monogenic OrphaCodes from TXT")

# Parse XML
tree = ET.parse(gene_assoc_xml)
root = tree.getroot()

# Containers
kept = []
excluded = []
kept_orphacodes = set()

# Define filter rules
valid_rules = {
    "36540": {"21436"},  # Group of disorders → Clinical group
    "36547": {"21394", "21422", "21415"},  # Disorder → Disease, Clinical syndrome, Morphological anomaly
    "36554": {"21450", "21443"}  # Subtype of disorder → Clinical subtype, Etiological subtype
}

# Iterate over disorders
for disorder in root.find(".//DisorderList"):
    orpha_code = disorder.findtext("OrphaCode", "").strip()
    if orpha_code not in monogenic_orphacodes:
        continue

    disorder_name = disorder.findtext("Name", "").strip()
    group_elem = disorder.find("DisorderGroup")
    type_elem = disorder.find("DisorderType")
    disorder_group_id = group_elem.attrib.get("id", "") if group_elem is not None else ""
    disorder_group_name = group_elem.findtext("Name", "").strip() if group_elem is not None else ""
    disorder_type_id = type_elem.attrib.get("id", "") if type_elem is not None else ""
    disorder_type_name = type_elem.findtext("Name", "").strip() if type_elem is not None else ""
    expert_link = disorder.findtext("ExpertLink", "").strip()

    # Association list count
    assoc_list = disorder.find("DisorderGeneAssociationList")
    assoc_count = assoc_list.attrib.get("count", "0") if assoc_list is not None else "0"

    # Iterate gene associations
    if assoc_list is not None:
        for gene_assoc in assoc_list.findall("DisorderGeneAssociation"):
            status = gene_assoc.findtext("DisorderGeneAssociationStatus/Name", "").strip()
            assoc_type = gene_assoc.findtext("DisorderGeneAssociationType/Name", "").strip()

            # Gene info
            gene = gene_assoc.find("Gene")
            gene_symbol = gene.findtext("Symbol", "").strip() if gene is not None else ""
            gene_name = gene.findtext("Name", "").strip() if gene is not None else ""
            gene_type = gene.findtext("GeneType/Name", "").strip() if gene is not None else ""
            gene_locus = ""
            if gene is not None:
                locus_list = gene.find("LocusList")
                if locus_list is not None and int(locus_list.attrib.get("count", "0")) > 0:
                    locus = locus_list.find("Locus/GeneLocus")
                    if locus is not None:
                        gene_locus = locus.text.strip()

            omim_ref = ""
            uniprot_ref = ""
            if gene is not None:
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
                orpha_code, disorder_name, disorder_group_name, disorder_group_id,
                disorder_type_name, disorder_type_id, gene_symbol, gene_name,
                gene_type, gene_locus, assoc_type, status,
                omim_ref, uniprot_ref, expert_link, assoc_count
            ]

            # Apply filter rules
            if disorder_group_id in valid_rules and disorder_type_id in valid_rules[disorder_group_id]:
                kept.append(row)
                kept_orphacodes.add(orpha_code)
            else:
                excluded.append(row)

# Debug prints
print(f"Total kept disorders: {len(kept_orphacodes)}")
print(f"Total excluded rows: {len(excluded)}")

# Save headers
headers = [
    "OrphaCode", "DisorderName", "DisorderGroup", "DisorderGroupID",
    "DisorderType", "DisorderTypeID", "GeneSymbol", "GeneName",
    "GeneType", "GeneLocus", "DisorderGeneAssociationType",
    "DisorderGeneAssociationStatus", "OMIM", "UniProt",
    "ExpertLink", "AssociationCount"
]

# Save kept CSV
kept_csv = os.path.join(output_folder, "kept_disorders.csv")
with open(kept_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(headers)
    writer.writerows(kept)
print(f"Saved kept_disorders.csv: {len(kept)} rows")

# Save excluded CSV
excluded_csv = os.path.join(output_folder, "excluded_disorders.csv")
with open(excluded_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(headers)
    writer.writerows(excluded)
print(f"Saved excluded_disorders.csv: {len(excluded)} rows")

# Save TXT of kept OrphaCodes
kept_txt = os.path.join(output_folder, "kept_orphacodes.txt")
with open(kept_txt, "w", encoding="utf-8") as f:
    for code in sorted(kept_orphacodes):
        f.write(f"{code}\n")
print(f"Saved TXT of kept OrphaCodes: {kept_txt}")

print("✅ Finished disorder classification pipeline")

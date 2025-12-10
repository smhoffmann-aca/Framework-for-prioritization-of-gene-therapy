import xml.etree.ElementTree as ET
import csv
import os

# File paths
tier1_matched_txt = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/tier1_matched_orphacodes.txt"
gene_assoc_xml = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Raw Datasets/genes_associated_ds.xml"
output_folder = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association"

os.makedirs(output_folder, exist_ok=True)

# Load Tier1 matched OrphaCodes
with open(tier1_matched_txt, "r") as f:
    tier1_matched_orphacodes = {line.strip() for line in f.readlines()}
print(f"Loaded {len(tier1_matched_orphacodes)} matched Tier1 OrphaCodes from TXT")

# Parse XML
tree = ET.parse(gene_assoc_xml)
root = tree.getroot()

# Store results
monogenic_rows = []
monogenic_orphacodes = set()

# Counters
count_monogenic = 0
count_polygenic = 0
count_no_assoc = 0

# Iterate disorders
for disorder in root.find(".//DisorderList"):
    orpha_code = disorder.findtext("OrphaCode", "").strip()
    if orpha_code not in tier1_matched_orphacodes:
        continue

    gene_assoc_list = disorder.find("DisorderGeneAssociationList")
    if gene_assoc_list is None:
        count_no_assoc += 1
        continue

    assoc_count = int(gene_assoc_list.attrib.get("count", "0"))
    if assoc_count == 1:
        count_monogenic += 1
    elif assoc_count > 1:
        count_polygenic += 1
    else:
        count_no_assoc += 1
        continue

    if assoc_count != 1:
        continue  # keep only monogenic disorders

    # Collect disorder-level info
    disorder_name = disorder.findtext("Name", "").strip()
    disorder_type_name = disorder.find("DisorderType/Name").text.strip() if disorder.find("DisorderType/Name") is not None else ""
    disorder_group_name = disorder.find("DisorderGroup/Name").text.strip() if disorder.find("DisorderGroup/Name") is not None else ""
    expert_link = disorder.findtext("ExpertLink", "").strip()

    # Collect gene association info
    for gene_assoc in gene_assoc_list.findall("DisorderGeneAssociation"):
        status = gene_assoc.findtext("DisorderGeneAssociationStatus/Name", "").strip()
        assoc_type = gene_assoc.findtext("DisorderGeneAssociationType/Name", "").strip()
        source_of_validation = gene_assoc.findtext("SourceOfValidation", "").strip()

        gene_symbol, gene_name, gene_type, gene_locus = "", "", "", ""
        omim_ref, uniprot_ref = "", ""

        gene = gene_assoc.find("Gene")
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
        monogenic_rows.append(row)
        monogenic_orphacodes.add(orpha_code)

# Debug prints
print("===== Debug Summary =====")
print(f"Total matched OrphaCodes: {len(tier1_matched_orphacodes)}")
print(f"Monogenic (AssocCount=1): {count_monogenic}")
print(f"Polygenic (AssocCount>1): {count_polygenic}")
print(f"No associations (AssocCount=0 or missing): {count_no_assoc}")

# Save TXT with OrphaCodes
monogenic_txt_file = os.path.join(output_folder, "monogenic_orphacodes.txt")
with open(monogenic_txt_file, "w", encoding="utf-8") as f:
    for code in sorted(monogenic_orphacodes):
        f.write(f"{code}\n")
print(f"Saved TXT of monogenic OrphaCodes: {monogenic_txt_file}")

# Save CSV with full info
csv_file = os.path.join(output_folder, "monogenic_associations.csv")
headers = [
    "OrphaCode", "DisorderName", "DisorderGroup", "DisorderType",
    "GeneSymbol", "GeneName", "GeneType", "GeneLocus",
    "DisorderGeneAssociationType", "DisorderGeneAssociationStatus",
    "SourceOfValidation", "OMIM", "UniProt", "ExpertLink", "AssocCount"
]

with open(csv_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(headers)
    writer.writerows(monogenic_rows)

print(f"Saved CSV of monogenic associations: {csv_file}")
print("✅ Finished monogenic association pipeline")

import xml.etree.ElementTree as ET
import csv
import os

# ===============================
# FILE PATHS
# ===============================

# Folder where the two Orphacode filter files are stored
filter_folder = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Gene Association/Monogenic Association/Disorder Classification/ClinGen/Natural History/Models/Clustering/Clinical Trials"

fileA = os.path.join(filter_folder, "tableA_orphacodes.txt")
fileB = os.path.join(filter_folder, "tableB_orphacodes.txt")

# Epidemiology XML dataset
xml_file = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders/Raw Datasets/epidemiology_ds.xml"

# Output CSV (saved in same folder as filter files)
csv_file = os.path.join(filter_folder, "epidemiology_filtered.csv")


# ===============================
# LOAD ORPHACODE FILTER LISTS
# ===============================
def load_orphacodes(path):
    with open(path, "r", encoding="utf-8") as f:
        return {line.strip() for line in f if line.strip()}

orphas_A = load_orphacodes(fileA)
orphas_B = load_orphacodes(fileB)

# Keep only codes present in BOTH lists
target_orphacodes = orphas_A.intersection(orphas_B)

print(f"Loaded {len(orphas_A)} Orphacodes from table A")
print(f"Loaded {len(orphas_B)} Orphacodes from table B")
print(f"→ {len(target_orphacodes)} Orphacodes present in BOTH lists\n")


# ===============================
# PARSE XML AND EXTRACT DATA
# ===============================
tree = ET.parse(xml_file)
root = tree.getroot()

# ===============================
# WRITE CSV
# ===============================
with open(csv_file, mode="w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file, delimiter=";")

    # Write header
    writer.writerow([
        "Orphacode", "Disorder Name", "Disorder Type", "Disorder Group",
        "Prevalence Source", "Prevalence Type", "Prevalence Qualification",
        "Prevalence Class", "Prevalence Value", "Geographic Area",
        "Validation Status", "Expert Link"
    ])

    # Iterate over disorder entries
    for disorder in root.findall(".//Disorder"):
        orpha_code = disorder.findtext("OrphaCode", "")

        # Only keep disorders that match BOTH filter lists
        if orpha_code not in target_orphacodes:
            continue

        disorder_name = disorder.findtext("Name", "")
        disorder_type = disorder.find("DisorderType/Name")
        disorder_type_name = disorder_type.text if disorder_type is not None else ""
        disorder_group = disorder.find("DisorderGroup/Name")
        disorder_group_name = disorder_group.text if disorder_group is not None else ""
        expert_link = disorder.findtext("ExpertLink", "")

        # Iterate through prevalence entries
        for prevalence in disorder.findall(".//Prevalence"):
            prevalence_source = prevalence.findtext("Source", "")
            prevalence_type = prevalence.find("PrevalenceType/Name")
            prevalence_type_name = prevalence_type.text if prevalence_type is not None else ""
            prevalence_qualification = prevalence.find("PrevalenceQualification/Name")
            prevalence_qualification_name = prevalence_qualification.text if prevalence_qualification is not None else ""
            prevalence_class = prevalence.findtext("PrevalenceClass/Name", "")
            prevalence_value = prevalence.findtext("ValMoy", "").replace(".", ",")
            prevalence_geographic = prevalence.find("PrevalenceGeographic/Name")
            prevalence_geographic_name = prevalence_geographic.text if prevalence_geographic is not None else ""
            prevalence_validation = prevalence.find("PrevalenceValidationStatus/Name")
            prevalence_validation_name = prevalence_validation.text if prevalence_validation is not None else ""

            # Write to CSV
            writer.writerow([
                orpha_code, disorder_name, disorder_type_name, disorder_group_name,
                prevalence_source, prevalence_type_name, prevalence_qualification_name,
                prevalence_class, prevalence_value, prevalence_geographic_name,
                prevalence_validation_name, expert_link
            ])

print(f"✅ CSV file created successfully:\n{csv_file}")

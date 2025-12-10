#!/usr/bin/env python3
import os
import csv
import xml.etree.ElementTree as ET

# === FILE PATHS ===
BASE_DIR = "/Users/sofiamorenohoffmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/M.Sc. Biomedical Sciences/Literature Review/Databases/Neurometabolic Disorders"

INPUT_XML = os.path.join(BASE_DIR, "Raw Datasets/natural_history_ds.xml")

FILTER_FILE = os.path.join(BASE_DIR, "Gene Association/Monogenic Association/Disorder Classification/ClinGen/definitive_orphacodes.txt")

OUTPUT_DIR = os.path.join(BASE_DIR, "Gene Association/Monogenic Association/Disorder Classification/ClinGen/Natural History")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === FILTER CRITERIA (IDs) ===
INCLUDE_IDS = {
    "23410",  # Autosomal dominant
    "23417",  # Autosomal recessive
    "23445",  # X-linked dominant
    "23431",  # X-linked recessive
    "23473",  # Y-linked
    "23466",  # Semi-dominant
}

EXCLUDE_IDS = {
    "23438",  # Mitochondrial inheritance
    "23424",  # Multigenic/multifactorial
    "23459",  # Oligogenic
    "23487",  # No data available
    "23494",  # Not applicable
    "23501",  # Not yet documented
}

# === FUNCTIONS ===
def load_filter_codes(filepath):
    """Load OrphaCodes from txt filter file."""
    codes = set()
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.isdigit():
                codes.add(line)
    return codes


def parse_inheritance(disorder):
    """Return list of inheritance IDs and names from a disorder node."""
    inheritance_nodes = disorder.findall(".//TypeOfInheritance")
    inheritance_ids = [n.attrib.get("id") for n in inheritance_nodes]
    inheritance_names = [n.findtext("Name[@lang='en']") for n in inheritance_nodes]
    return inheritance_ids, inheritance_names


def should_include(inheritance_ids):
    """Inclusion rule: keep if has at least one include ID and not all excluded."""
    has_include = any(i in INCLUDE_IDS for i in inheritance_ids)
    all_excluded = all(i in EXCLUDE_IDS for i in inheritance_ids) if inheritance_ids else True
    return has_include and not all_excluded


def parse_dataset(xml_file, filter_codes):
    """Parse XML and return included and excluded disorder dicts."""
    tree = ET.parse(xml_file)
    root = tree.getroot()

    all_disorders = root.findall(".//Disorder")
    total_disorders = len(all_disorders)

    included, excluded = [], []
    matched_codes = set()

    for disorder in all_disorders:
        orphacode = disorder.findtext("OrphaCode")
        if not orphacode:
            continue

        if orphacode in filter_codes:
            matched_codes.add(orphacode)

            name = disorder.findtext("Name[@lang='en']")
            dgroup = disorder.findtext("DisorderGroup/Name[@lang='en']")
            dtype = disorder.findtext("DisorderType/Name[@lang='en']")

            onset_list = [n.findtext("Name[@lang='en']") for n in disorder.findall(".//AverageAgeOfOnset")]
            onset = ";".join(onset_list) if onset_list else "No data"

            inh_ids, inh_names = parse_inheritance(disorder)
            inh_count = len(inh_ids)
            inh_names_str = ";".join(inh_names) if inh_names else "No data"
            inh_ids_str = ";".join(inh_ids) if inh_ids else "No data"

            row = {
                "OrphaCode": orphacode,
                "Name": name,
                "DisorderGroup": dgroup,
                "DisorderType": dtype,
                "AverageAgeOfOnset": onset,
                "InheritanceCount": inh_count,
                "InheritanceTypes": inh_names_str,
                "InheritanceIDs": inh_ids_str,
            }

            if should_include(inh_ids):
                included.append(row)
            else:
                excluded.append(row)

    return total_disorders, included, excluded, matched_codes


def save_csv(rows, filepath):
    """Save disorder rows to CSV."""
    if not rows:
        return
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def save_txt(codes, filepath):
    """Save OrphaCodes to TXT."""
    with open(filepath, "w", encoding="utf-8") as f:
        for c in sorted(codes, key=lambda x: int(x)):
            f.write(c + "\n")


# === MAIN ===
def main():
    filter_codes = load_filter_codes(FILTER_FILE)
    print(f"🧩 Loaded {len(filter_codes)} OrphaCodes from filter file: {FILTER_FILE}")

    total_disorders, included, excluded, matched_codes = parse_dataset(INPUT_XML, filter_codes)

    print("\n=== SUMMARY ===")
    print(f"Total disorders in natural history dataset: {total_disorders}")
    print(f"Filter OrphaCodes provided: {len(filter_codes)}")
    print(f"OrphaCodes found in dataset: {len(matched_codes)}")
    print(f"Included after inheritance filtering: {len(included)}")
    print(f"Excluded after inheritance filtering: {len(excluded)}")

    # Save CSVs
    save_csv(included, os.path.join(OUTPUT_DIR, "included_inheritance_disorders.csv"))
    save_csv(excluded, os.path.join(OUTPUT_DIR, "excluded_inheritance_disorders.csv"))

    # Save TXT lists
    save_txt([r["OrphaCode"] for r in included], os.path.join(OUTPUT_DIR, "included_inheritance_orphacodes.txt"))
    save_txt([r["OrphaCode"] for r in excluded], os.path.join(OUTPUT_DIR, "excluded_inheritance_orphacodes.txt"))

    print(f"\n✅ Output files saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

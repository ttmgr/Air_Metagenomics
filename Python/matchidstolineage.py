import glob
import os

def load_lineage_from_file(lineage_file):
    lineage_dict = {}
    with open(lineage_file, "r") as infile:
        for line in infile:
            fields = line.strip().split("\t")
            taxid = fields[2]
            lineage = fields[3]
            lineage_dict[taxid] = lineage
    print(f"Loaded {len(lineage_dict)} lineages.")
    return lineage_dict

def process_file(input_file, lineage_dict):
    basename = os.path.splitext(input_file)[0]
    output_file = f"{basename}_with_lineage_added.txt"
    
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        lines_processed = 0
        for line in infile:
            fields = line.strip().split("\t")
            taxid = fields[-1]
            lineage = lineage_dict.get(taxid, "Lineage not found")
            if lineage != "Lineage not found":
                outfile.write(f"{line.strip()}\t{lineage}\n")
                lines_processed += 1
        print(f"Processed {lines_processed} lines in {input_file}")

if __name__ == "__main__":
    lineage_dict = load_lineage_from_file("output_with_full_lineage.txt")
    
    for dmnd_file in glob.glob("*.dmnd_out"):
        process_file(dmnd_file, lineage_dict)

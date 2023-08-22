import csv
import glob

def process_file(file_path, database):
    # Read the file and extract the required columns
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            subject_seq_id = row['subject_seq_id']
            organism = row['organism']
            taxonomy = row['taxonomy']

            # Add to the database if not already present
            if subject_seq_id not in database:
                database[subject_seq_id] = {'organism': organism, 'taxonomy': taxonomy}

# Initialize an empty database
database = {}

# Process all .tsv files in the current directory
for file_path in glob.glob('*.tsv'):
    process_file(file_path, database)

# Write the database to a new TSV file
output_file_path = 'database.tsv'
with open(output_file_path, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['subject_seq_id', 'organism', 'taxonomy'])
    for subject_seq_id, data in database.items():
        writer.writerow([subject_seq_id, data['organism'], data['taxonomy']])

print(f"Database written to {output_file_path}")

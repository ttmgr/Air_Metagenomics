# Define the absolute path to the configuration file
CONFIG_FILE="/home/haicu/ttreska57/dorado/dorado-0.5.0-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v4.3.0"

# Step 1: Run dorado basecaller
echo "Running dorado basecaller..."
${DORADO_BIN} duplex  ${CONFIG_FILE} -r /lustre/groups/hpc/urban_lab/backup/barcelona_air/valle_dhebron/ > basecalled_duplex.bam -t 2

# Check if basecalling was successful
if [ $? -ne 0 ]; then
    echo "Error in basecalling step."
    exit 1
fi

echo "Process completed successfully."

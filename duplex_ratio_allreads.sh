# Define your input file and output file
input_file="./basecalled_duplex.bam"
output_file="./read_ratio_results.txt"

# Count simplex reads (tags: dx:i:0 or dx:i:-1)
simplex_count=$(samtools view $input_file | grep -E 'dx:i:0|dx:i:-1' | wc -l)

# Count duplex reads (tag: dx:i:1)
duplex_count=$(samtools view $input_file | grep 'dx:i:1' | wc -l)

# Calculate the ratio of duplex to simplex reads
if [ $duplex_count -gt 0 ]; then
  ratio=$(echo "scale=2; $duplex_count / $simplex_count" | bc)
  else
    ratio="No duplex reads found."
    fi

    # Save the results to a file
    echo "Simplex reads: $simplex_count" > $output_file
    echo "Duplex reads: $duplex_count" >> $output_file
    echo "Ratio of Simplex to Duplex Reads: $ratio" >> $output_file

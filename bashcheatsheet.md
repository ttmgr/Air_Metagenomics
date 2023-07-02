```shell

# FILE AND DIRECTORY MANAGEMENT

# List files and directories in the current directory
ls

# List files and directories in the current directory with details
ls -l

# List all files and directories (including hidden ones) in the current directory
ls -a

# Change directory to a specified path
cd /path/to/directory

# Go back to the previous directory
cd -

# Create a new directory
mkdir new_directory_name

# Remove an empty directory
rmdir directory_name

# Remove a non-empty directory and its contents
rm -r directory_name

# Remove a file
rm file_name

# Copy a file or directory
cp source destination

# Move (or rename) a file or directory
mv source destination

# Create a symbolic link
ln -s target link_name

# FILE CONTENT MANIPULATION

# Print the contents of a file
cat file_name

# Print the first 10 lines of a file
head file_name

# Print the last 10 lines of a file
tail file_name

# Search for a pattern in a file
grep 'pattern' file_name

# Search and replace a pattern in a file
sed 's/old_pattern/new_pattern/g' input_file > output_file

# PROCESS MANAGEMENT

# List currently running processes
ps

# List all currently running processes with details
ps -aux

# Stop a process by process ID (PID)
kill PID

# Stop a process by name
killall process_name

# SYSTEM INFORMATION

# Show system information
uname -a

# Show disk usage
df -h

# Show memory usage
free -h

# Show system uptime and load
uptime

# NETWORKING

# Display network interfaces and their configurations
ifconfig

# Display routing table
route

# Test network connectivity
ping host_name_or_ip_address

# Download a file from the internet
wget URL

# Transfer files over the network using the SCP protocol
scp source user@host:destination

# PERMISSIONS AND OWNERSHIP

# Change file permissions (chmod)
chmod mode file_name

# Change file owner and group (chown)
chown owner:group file_name

# ARCHIVE AND COMPRESS

# Create a tarball archive
tar -cvf archive_name.tar file_or_directory

# Extract a tarball archive
tar -xvf archive_name.tar

# Create a compressed gzip archive
tar -cvzf archive_name.tar.gz file_or_directory

# Extract a compressed gzip archive
tar -xvzf archive_name.tar.gz

# FILE SORTING AND COMPARISON

# Sort a file by lines
sort input_file > sorted_file

# Sort a file numerically by the first column
sort -n -k1 input_file > sorted_file

# Sort a file in reverse order
sort -r input_file > sorted_file

# Remove duplicate lines from a sorted file
uniq sorted_input_file > unique_file

# Compare two files line by line
diff file1 file2

# FILE SLICING AND COLUMN EXTRACTION

# Print specific lines from a file
sed -n 'start_line,end_line p' input_file > output_file

# Remove specific lines from a file
sed 'start_line,end_line d' input_file > output_file

# Extract specific columns from a file (e.g., columns 1 and 3)
cut -f1,3 input_file > output_file

# Combine multiple files by columns
paste file1 file2 > combined_file

# DATA PROCESSING AND FILTERING

# Summarize data (count, sum, mean, min, max)
awk '{sum+=$1; count+=1} END {print "Count:", count, "Sum:", sum, "Mean:", sum/count, "Min:", min, "Max:", max}' input_file

# Filter lines based on a condition (e.g., first column greater than 10)
awk '$1 > 10' input_file > output_file

# Filter lines containing a specific pattern
grep 'pattern' input_file > output_file

# Filter lines not containing a specific pattern
grep -v 'pattern' input_file > output_file

# Merge lines based on a common field (e.g., merge on the first field)
join -1 1 -2 1 file1 file2 > output_file

# Calculate the number of lines, words, and characters in a file
wc input_file

# Calculate the number of lines in a file
wc -l input_file

# Calculate the number of words in a file
wc -w input_file

# Calculate the number of characters in a file
wc -c input_file

# FINDING AND SEARCHING

# Find files in a directory that match a pattern
find /path/to/search -name 'pattern'

# Find files in a directory that are larger than a specific size (e.g., 100MB)
find /path/to/search -type f -size +100M

# Replace a pattern in a file
sed 's/old_pattern/new_pattern/g' input_file > output_file

# Replace a pattern in a file in-place
sed -i 's/old_pattern/new_pattern/g' input_file

# PERL ONE-LINERS

# Sum a column of numbers
perl -lane '$sum += $F[0]; END {print $sum}' input_file

# Find the median of a column of numbers
perl -lane 'push @array, $F[0]; END { @sorted = sort { $a <=> $b } @array; $median = @sorted % 2 ? $sorted[(@sorted - 1) / 2] : ($sorted[@sorted / 2 - 1] + $sorted[@sorted / 2]) / 2; print $median; }' input_file

# Remove duplicate lines while preserving the original order
perl -ne 'print if !$seen{$_}++' input_file > unique_file

# TEXT AND DATA TRANSFORMATION

# Transpose rows to columns and columns to rows
awk '{ for (i=1; i<=NF; i++) a[NR][i] = $i; max_nf = (NF > max_nf ? NF : max_nf); } END { for (j=1; j<=max_nf; j++) { str=a[1][j]; for (i=2; i<=NR; i++) str=str" "a[i][j]; print str; } }' input_file > output_file

# Convert CSV file to TSV (Tab-Separated Values)
awk -F, '{gsub(/,/,"\t"); print}' input_file > output_file

# Convert TSV file to CSV (Comma-Separated Values)
awk -F'\t' '{gsub(/\t/,","); print}' input_file > output_file

# PARALLEL PROCESSING

# Run multiple commands in parallel
parallel ::: 'command1' 'command2' 'command3'

# Run a command for each file in a directory in parallel
ls input_directory | parallel 'command {} > output_directory/{}.out'

# Download multiple files in parallel
cat urls.txt | parallel -j4 wget -P /path/to/save

# Calculate sum of a column in a file
awk '{sum += $1} END {print sum}' input_file.txt

# Calculate mean (average) of a column in a file
awk '{sum += $1; n++} END {if (n > 0) print sum / n; else print "No data"}' input_file.txt

# Calculate median of a column in a file
sort -n input_file.txt | awk 'BEGIN {c = 0} {a[c++] = $1} END {if (c % 2 == 1) print a[int(c / 2)]; else print (a[c / 2] + a[c / 2 - 1]) / 2}'

# Calculate mode of a column in a file
awk '{count[$1]++} END {for (i in count) if (count[i] > max) {max = count[i]; mode = i}} END {print mode}' input_file.txt

# Calculate variance of a column in a file
awk '{sum += $1; sumsq += $1*$1} END {print (sumsq - (sum*sum) / NR) / NR}' input_file.txt

# Calculate standard deviation of a column in a file
awk '{sum += $1; sumsq += $1*$1} END {print sqrt((sumsq - (sum*sum) / NR) / NR)}' input_file.txt




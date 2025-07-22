# Troubleshooting Guide

This guide provides solutions to common issues you might encounter while setting up or running the nanopore metagenomics pipeline.

---

## ⚠️ Issue: "Command not found" or Tool Incompatibility Errors

This is the most common issue and usually relates to the software environment.

* **Symptom**: The pipeline fails immediately with an error like `bash_scripts/01_read_processing.sh: line 25: porechop: command not found`.
* **Cause**: The Conda environment (`nanopore-metagenomics`) is not activated, or it was not created correctly.
* **Solution**:
    1.  **Activate the Environment**: Ensure you have activated the correct environment before running the pipeline.
        ```bash
        mamba activate nanopore-metagenomics
        ```
    2.  **Verify Installation**: If the error persists, the environment may be corrupted or incomplete. The best solution is to rebuild it. First, remove the old environment, then create it again using the `environment.yaml` file.
        ```bash
        # Deactivate if currently active
        mamba deactivate

        # Remove the broken environment
        mamba env remove --name nanopore-metagenomics

        # Re-create it from the YAML file
        mamba env create -f env/environment.yaml
        ```

---

## 🗃️ Issue: Database Errors

Problems with database paths or files are another frequent cause of pipeline failure.

* **Symptom**: A tool (like Kraken2, Bakta, or AMRFinderPlus) fails with an error message like `database not found`, `Cannot open file`, or `Error reading database`.
* **Cause**: The paths configured in `run_pipeline.sh` are incorrect, or the database files are corrupted or have incorrect permissions.
* **Solutions**:
    1.  **Check Paths in `run_pipeline.sh`**: This is the most likely cause. Open `bash_scripts/run_pipeline.sh` and verify that the database paths (`KRAKEN2_DB_PATH`, `BAKTA_DB_PATH`, etc.) **exactly match** the locations where you downloaded the databases with `download_databases.sh`.
    2.  **File Permissions**: Make sure you have read permissions for the database files. You can check this with `ls -l /path/to/your/databases`.
    3.  **Use Pre-built Kraken2 Databases**: The `download_databases.sh` script builds the Kraken2 database from source, which is very slow and memory-intensive. A faster and more reliable alternative is to download a pre-built database index.
        * You can find a list of available pre-built indexes here: [https://benlangmead.github.io/aws-indexes/](https://benlangmead.github.io/aws-indexes/)
        * **Example**: To download and use the 8GB Standard database:
            ```bash
            # Navigate to your main database directory
            cd /path/to/your/databases

            # Download the pre-built index
            wget [https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20230314.tar.gz](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20230314.tar.gz)

            # Extract the archive
            tar -zxvf k2_standard_8gb_20230314.tar.gz

            # Clean up the archive file
            rm k2_standard_8gb_20230314.tar.gz
            ```
        * After extracting, update the `KRAKEN2_DB_PATH` in `run_pipeline.sh` to point to the new directory (e.g., `/path/to/your/databases/k2_standard_8gb_20230314`).

---

## 💥 Issue: Pipeline Fails on a Specific Step

* **Symptom**: The pipeline runs for a while but then stops unexpectedly during one of the stages (e.g., during assembly with Flye).
* **Cause**: This could be due to resource limitations (memory/CPU), an issue with a specific input file, or a tool-specific bug.
* **Solutions**:
    1.  **Check the Log Files**: The main output directory contains detailed logs. The main log is `pipeline.log`, but each stage may also have its own log file (e.g., `assembly.log`). Open the relevant log file and scroll to the bottom to find the specific error message.
    2.  **Check Input Files**: Ensure your input FASTQ files in the `INPUT_FASTQ_DIR` are not empty or corrupted. A common mistake is an incorrect path, leading to the script finding no files to process.
    3.  **Resource Issues**: Metagenome assembly (Flye) and binning (MetaWRAP) are very memory-intensive. If the pipeline crashes without a clear error in the log, it might be an "Out of Memory" (OOM) error from the system.
        * Try running the pipeline on a machine with more RAM.
        * Reduce the number of `THREADS` in `run_pipeline.sh`.
        * On a Linux system, you can check for OOM errors using the command `dmesg -T | grep -i "out of memory"`.

---

## 📉 Issue: Low Read Counts or Poor Quality Results

* **Symptom**: The pipeline finishes, but the results are poor. For example, NanoStat shows very few reads passed filtering, or Kraken2 classifies almost nothing.
* **Cause**: This usually points to issues with the initial data quality or the filtering parameters.
* **Solutions**:
    1.  **Review Pre-QC Data**: Go back to your basecalling and demultiplexing results. A high number of reads in the `unclassified` folder from your basecaller indicates a problem with the barcoding or demultiplexing step.
    2.  **Check for Contamination**: If a large percentage of your reads classify as *Homo sapiens* or other known contaminants, you may need to add a pre-filtering step to remove these reads before starting the main pipeline.
    3.  **Adjust Filtering Parameters**: The `01_read_processing.sh` script uses default filtering parameters (`-q 9` and `-l 500`). These are a good starting point but may not be optimal for every dataset. If you are losing too many reads, consider lowering the quality threshold (e.g., `-q 8`) or the minimum length (e.g., `-l 300`). Conversely, if your downstream results are poor, you may need to increase the stringency.

If you encounter tool-specific errors not covered here, the best course of action is to copy the error message and consult the official documentation or GitHub page for that particular tool.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consolidated Python script for Nanopore Metagenomics Analysis.

This script combines functionalities for:
1.  Extracting read quality metrics from NanoStat.
2.  Processing and summarizing Kraken2 taxonomic classification reports.
3.  Visualizing data, including read length histograms and read count plots.
4.  Generating a comprehensive HTML report.
"""
import os
import pandas as pd
import argparse
import logging
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from jinja2 import Template
from collections import defaultdict

# --- Basic Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Nanostat Metrics Extraction ---
def extract_nanostat_metrics(directory):
    """Extracts metrics from NanoStat report files."""
    data_list = []
    logger.info(f"Scanning directory for NanoStat files: {directory}")
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            try:
                barcode = filename.split('_')[1] # Assumes format 'filtered_barcode01_... .txt'
                file_stats = {'Barcode': barcode, 'Mean_Read_Length': np.nan, 'Median_Read_Length': np.nan, 'Read_Length_N50': np.nan, 'Number_of_Reads': np.nan}
                with open(os.path.join(directory, filename), 'r') as file:
                    for line in file:
                        if 'Mean read length:' in line:
                            file_stats['Mean_Read_Length'] = float(line.split(':')[1].strip().replace(',', ''))
                        elif 'Median read length:' in line:
                            file_stats['Median_Read_Length'] = float(line.split(':')[1].strip().replace(',', ''))
                        elif 'Read length N50:' in line:
                            file_stats['Read_Length_N50'] = float(line.split(':')[1].strip().replace(',', ''))
                        elif 'Number of reads:' in line:
                            file_stats['Number_of_Reads'] = float(line.split(':')[1].strip().replace(',', ''))
                data_list.append(file_stats)
            except IndexError:
                logger.warning(f"Could not extract barcode from filename: {filename}. Skipping.")
                continue
    if not data_list:
        logger.warning("No NanoStat files were processed. Please check the directory and file naming.")
    return pd.DataFrame(data_list)

# --- Read Length Histogram ---
def plot_read_length_histogram(fastq_path, output_path):
    """Generates and saves a read length histogram from a FASTQ file."""
    logger.info(f"Extracting read lengths from {fastq_path}")
    read_lengths = []
    with open(fastq_path, 'r') as file:
        for i, line in enumerate(file):
            if i % 4 == 1:  # Sequence line
                read_lengths.append(len(line.strip()))
    
    if not read_lengths:
        logger.warning(f"No reads found in {fastq_path}. Cannot generate histogram.")
        return

    read_lengths_array = np.array(read_lengths)
    
    plt.figure(figsize=(12, 7))
    plt.hist(read_lengths_array, bins=np.logspace(np.log10(max(1, np.min(read_lengths_array))), np.log10(np.max(read_lengths_array)), 80), color='lightgreen', edgecolor='black')
    plt.xscale('log')
    
    median_read_length = np.median(read_lengths_array)
    plt.axvline(median_read_length, color='black', linestyle='--', label=f'Median: {median_read_length:.2f} bp')
    
    plt.xlabel('Read Length (bp, log scale)')
    plt.ylabel('Number of Reads')
    plt.title(f'Read Length Distribution for {os.path.basename(fastq_path)}')
    plt.legend()
    plt.grid(axis='y', alpha=0.75)
    
    output_filename = os.path.join(output_path, f"{Path(fastq_path).stem}_length_histogram.png")
    plt.savefig(output_filename)
    plt.close()
    logger.info(f"Read length histogram saved to {output_filename}")

# --- Kraken2 Report Processing ---
def process_single_kraken_report(file_path):
    """Processes a single Kraken2 report to calculate relative abundances."""
    df = pd.read_csv(file_path, sep='\t', header=None, names=["Percentage", "Reads", "Reads_at_Level", "Rank", "TaxID", "Name"])
    df['Name'] = df['Name'].str.strip()
    
    # Filter out unclassified and human reads
    df_filtered = df[(df['Rank'] != 'U') & (df['Name'] != 'Homo')]
    df_taxa = df_filtered[df_filtered['Rank'].isin(['P', 'G'])] # Phylum and Genus
    
    # Calculate relative abundance for P and G ranks
    taxa_reads = df_taxa.groupby(['Rank', 'Name'])['Reads'].sum()
    total_reads_by_rank = taxa_reads.groupby('Rank').sum()
    
    if total_reads_by_rank.empty:
        return pd.DataFrame(columns=['Rank.code', 'Name', 'Relative_Abundance'])

    relative_abundances = taxa_reads.div(total_reads_by_rank, level='Rank').fillna(0) * 100
    df_relative = relative_abundances.reset_index()
    df_relative.columns = ['Rank.code', 'Name', 'Relative_Abundance']
    
    return df_relative[df_relative['Relative_Abundance'] >= 1.0]

def process_all_kraken_reports(directory_path):
    """Processes all Kraken2 reports in a directory for relative abundance."""
    combined_df = pd.DataFrame()
    logger.info(f"Scanning directory for Kraken2 reports: {directory_path}")
    for filename in os.listdir(directory_path):
        if filename.endswith(".txt"):
            file_path = os.path.join(directory_path, filename)
            df = process_single_kraken_report(file_path)
            df['Sample'] = Path(filename).stem.replace('report_', '')
            combined_df = pd.concat([combined_df, df], ignore_index=True)
    if combined_df.empty:
        logger.warning("No Kraken2 reports were processed. Please check the directory.")
    return combined_df

def extract_kraken_read_counts(directory_path):
    """Extracts total phylum and genus read counts from Kraken2 reports."""
    data_list = []
    logger.info(f"Scanning directory for Kraken2 reports to count reads: {directory_path}")
    for filename in os.listdir(directory_path):
        if filename.endswith(".txt"):
            try:
                # Assuming format 'report_filtered_barcode01_passed.txt'
                barcode = filename.split('_')[2]
                phylum_sum = 0
                genus_sum = 0
                with open(os.path.join(directory_path, filename), 'r') as file:
                    for line in file:
                        parts = line.split('\t')
                        if len(parts) >= 6:
                            rank = parts[3].strip()
                            count = int(parts[1].strip())
                            if rank == 'P':
                                phylum_sum += count
                            elif rank == 'G':
                                genus_sum += count
                data_list.append({'Barcode': barcode, 'Phylum_Reads': phylum_sum, 'Genus_Reads': genus_sum})
            except IndexError:
                logger.warning(f"Could not extract barcode from filename: {filename} using the expected format. Skipping.")
                continue

    if not data_list:
        logger.warning("No Kraken2 reports were processed for read counts.")
    return pd.DataFrame(data_list)
    
# --- Plotting Classified Reads ---
def plot_classified_read_counts(csv_path, output_dir):
    """Creates bar plots of classified read counts from a summary CSV file."""
    try:
        data = pd.read_csv(csv_path)
        # Infer location from sample ID
        data['Location'] = data['Sample ID'].str.extract(r'(^.*?)(?:_\d| \d)')[0].fillna('Unknown')
    except Exception as e:
        logger.error(f"Failed to load or process CSV {csv_path}: {e}")
        return

    sns.set(style="whitegrid")
    for location in data['Location'].unique():
        location_data = data[data['Location'] == location].copy()
        # Create a shorter sample ID for plotting
        location_data['Short Sample ID'] = location_data['Sample ID'].str.replace(f'{location}', '', n=1).str.strip(' _')
        
        plt.figure(figsize=(14, 10))
        
        plot_cols = [col for col in ['Reads after filtering', 'Reads mapped to phylum', 'Reads mapped to genus'] if col in location_data.columns]
        if not plot_cols:
            logger.warning(f"No data columns to plot for location: {location}")
            continue

        ax = location_data.set_index('Short Sample ID')[plot_cols].plot(kind='bar')
        
        ax.set_title(f'Read Counts for {location}')
        ax.set_xlabel('Sample Replicate')
        ax.set_ylabel('Number of Reads')
        
        if 'Reads after filtering' in location_data.columns:
            ylim_max = location_data['Reads after filtering'].max() * 1.1
            ax.set_ylim(0, ylim_max)
        
        ax.legend(title='Metrics')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        output_filename = os.path.join(output_dir, f"{location}_read_counts.png")
        plt.savefig(output_filename)
        plt.close()
        logger.info(f"Read count plot saved to {output_filename}")

# --- HTML Report Generation ---
class MetagenomicsReportGenerator:
    """Generate comprehensive HTML reports for metagenomics analysis."""
    
    def __init__(self, input_dir, output_file):
        self.input_dir = Path(input_dir)
        self.output_file = Path(output_file)
        self.data = defaultdict(dict)
        self.figures = {}
        logger.info(f"ReportGenerator initialized for project: {self.input_dir}")

    def run(self):
        """Execute all steps to generate the report."""
        self._collect_data()
        self._create_visualizations()
        self._generate_report()

    def _collect_data(self):
        """Collect all necessary data from the project directory."""
        logger.info("Collecting data for report...")
        # Placeholder for collecting various data sources
        # For example, collecting nanostat summary:
        nanostat_dir = self.input_dir / "processing" / "nanostat"
        if nanostat_dir.exists():
            self.data['preprocessing_stats'] = extract_nanostat_metrics(nanostat_dir)
        
        kraken_dir = self.input_dir / "processing" / "kraken2_read_classification"
        if kraken_dir.exists():
            self.data['kraken_abundance'] = process_all_kraken_reports(kraken_dir)

    def _create_visualizations(self):
        """Create Plotly visualizations."""
        logger.info("Creating visualizations...")
        if not self.data['preprocessing_stats'].empty:
            self._plot_preprocessing_stats()
        if 'kraken_abundance' in self.data and not self.data['kraken_abundance'].empty:
            self._plot_taxonomic_heatmap()

    def _plot_preprocessing_stats(self):
        df = self.data['preprocessing_stats']
        fig = make_subplots(rows=2, cols=2, subplot_titles=('Read Counts', 'Total Bases', 'Mean Read Length', 'Read Length N50'))
        fig.add_trace(go.Bar(x=df['Barcode'], y=df['Number_of_Reads'], name='Reads'), row=1, col=1)
        fig.add_trace(go.Bar(x=df['Barcode'], y=df['Number_of_Reads']*df['Mean_Read_Length'], name='Bases'), row=1, col=2)
        fig.add_trace(go.Bar(x=df['Barcode'], y=df['Mean_Read_Length'], name='Mean Length'), row=2, col=1)
        fig.add_trace(go.Bar(x=df['Barcode'], y=df['Read_Length_N50'], name='N50'), row=2, col=2)
        fig.update_layout(height=800, title_text="Read Quality Metrics Overview", showlegend=False)
        self.figures['preprocessing_plot'] = fig.to_html(full_html=False, include_plotlyjs='cdn')

    def _plot_taxonomic_heatmap(self):
        df = self.data['kraken_abundance']
        top_taxa = df.groupby('Name')['Relative_Abundance'].sum().nlargest(25).index
        df_top = df[df['Name'].isin(top_taxa)]
        heatmap_data = df_top.pivot_table(index='Name', columns='Sample', values='Relative_Abundance').fillna(0)
        fig = px.imshow(heatmap_data, labels=dict(x="Sample", y="Taxon", color="Abundance (%)"), title="Top 25 Most Abundant Taxa (Genus/Phylum)")
        self.figures['taxonomic_heatmap'] = fig.to_html(full_html=False, include_plotlyjs='cdn')

    def _generate_report(self):
        """Render the Jinja2 template with collected data and plots."""
        template_str = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Metagenomics Analysis Report</title>
            <style>
                body { font-family: sans-serif; margin: 2em; }
                h1, h2 { color: #333; }
                .content { max-width: 1200px; margin: auto; }
                .plot { margin-top: 2em; }
            </style>
        </head>
        <body>
            <div class="content">
                <h1>Metagenomics Pipeline Report</h1>
                <p><strong>Report generated on:</strong> {{ generation_date }}</p>
                <p><strong>Project Directory:</strong> {{ project_dir }}</p>
                
                <h2>Preprocessing Metrics</h2>
                {% if preprocessing_table %}
                    {{ preprocessing_table }}
                {% else %}
                    <p>No preprocessing data found.</p>
                {% endif %}
                <div class="plot">{{ figures.get('preprocessing_plot', '') }}</div>
                
                <h2>Taxonomic Classification</h2>
                {% if kraken_table %}
                    {{ kraken_table }}
                {% else %}
                     <p>No Kraken2 abundance data found.</p>
                {% endif %}
                <div class="plot">{{ figures.get('taxonomic_heatmap', '') }}</div>
            </div>
        </body>
        </html>
        """
        template = Template(template_str)

        # Prepare tables for rendering
        preprocessing_table_html = self.data['preprocessing_stats'].to_html(index=False, classes='table') if not self.data['preprocessing_stats'].empty else None
        kraken_table_html = self.data['kraken_abundance'].to_html(index=False, classes='table') if 'kraken_abundance' in self.data and not self.data['kraken_abundance'].empty else None

        html_content = template.render(
            generation_date=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            project_dir=self.input_dir,
            figures=self.figures,
            preprocessing_table=preprocessing_table_html,
            kraken_table=kraken_table_html
        )
        
        with open(self.output_file, 'w') as f:
            f.write(html_content)
        logger.info(f"HTML report saved to {self.output_file}")


# --- Main Function and Argument Parsing ---
def main():
    parser = argparse.ArgumentParser(
        description="A consolidated script for nanopore metagenomics data analysis.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='task', required=True, help='Available tasks')

    # Sub-parser for nanostat
    p_ns = subparsers.add_parser('nanostat', help='Extract metrics from NanoStat output files.')
    p_ns.add_argument('--input_dir', required=True, help='Directory containing NanoStat .txt files.')
    p_ns.add_argument('--output_file', required=True, help='Path to save the output CSV file.')

    # Sub-parser for read length histogram
    p_hist = subparsers.add_parser('length_histogram', help='Generate a read length histogram from a FASTQ file.')
    p_hist.add_argument('--input_file', required=True, help='Path to the input FASTQ file.')
    p_hist.add_argument('--output_dir', required=True, help='Directory to save the output histogram PNG.')
    
    # Sub-parser for Kraken2 relative abundance
    p_kr = subparsers.add_parser('kraken_abundance', help='Process Kraken2 reports for relative abundance of major taxa.')
    p_kr.add_argument('--input_dir', required=True, help='Directory containing Kraken2 report .txt files.')
    p_kr.add_argument('--output_file', required=True, help='Path to save the output CSV file.')

    # Sub-parser for Kraken2 read counts
    p_kc = subparsers.add_parser('kraken_counts', help='Extract total phylum and genus read counts from Kraken2 reports.')
    p_kc.add_argument('--input_dir', required=True, help='Directory containing Kraken2 report .txt files.')
    p_kc.add_argument('--output_file', required=True, help='Path to save the output CSV file.')
    
    # Sub-parser for plotting read counts
    p_plot = subparsers.add_parser('plot_counts', help='Plot classified read counts from a summary CSV file.')
    p_plot.add_argument('--input_file', required=True, help='Path to the summary CSV file.')
    p_plot.add_argument('--output_dir', required=True, help='Directory to save the output plots.')
    
    # Sub-parser for generating the final report
    p_report = subparsers.add_parser('create_report', help='Generate a final HTML report from pipeline results.')
    p_report.add_argument('--input_dir', required=True, help='Main pipeline output directory containing all results.')
    p_report.add_argument('--output_file', required=True, help='Path for the final HTML report file.')

    args = parser.parse_args()

    # --- Task Execution ---
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True) if hasattr(args, 'output_file') else os.makedirs(args.output_dir, exist_ok=True)

    if args.task == 'nanostat':
        df = extract_nanostat_metrics(args.input_dir)
        df.to_csv(args.output_file, index=False)
        logger.info(f"NanoStat metrics summary saved to {args.output_file}")
    
    elif args.task == 'length_histogram':
        plot_read_length_histogram(args.input_file, args.output_dir)

    elif args.task == 'kraken_abundance':
        df = process_all_kraken_reports(args.input_dir)
        df.to_csv(args.output_file, index=False)
        logger.info(f"Kraken2 relative abundance summary saved to {args.output_file}")
        
    elif args.task == 'kraken_counts':
        df = extract_kraken_read_counts(args.input_dir)
        df.to_csv(args.output_file, index=False)
        logger.info(f"Kraken2 read count summary saved to {args.output_file}")

    elif args.task == 'plot_counts':
        plot_classified_read_counts(args.input_file, args.output_dir)
    
    elif args.task == 'create_report':
        report_generator = MetagenomicsReportGenerator(args.input_dir, args.output_file)
        report_generator.run()

if __name__ == '__main__':
    main()

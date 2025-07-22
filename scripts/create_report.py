#!/usr/bin/env python3
"""
Comprehensive Report Generation for Nanopore Metagenomics Pipeline
Generates HTML reports with interactive visualizations
"""

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from jinja2 import Template
import argparse
import logging
from collections import defaultdict

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MetagenomicsReportGenerator:
    """Generate comprehensive HTML reports for metagenomics analysis"""
    
    def __init__(self, input_dir, output_file):
        self.input_dir = Path(input_dir)
        self.output_file = Path(output_file)
        self.data = defaultdict(dict)
        self.figures = {}
        
    def collect_data(self):
        """Collect all analysis data from the pipeline output"""
        logger.info("Collecting data from pipeline output...")
        
        # Collect preprocessing stats
        self._collect_preprocessing_stats()
        
        # Collect assembly stats
        self._collect_assembly_stats()
        
        # Collect classification results
        self._collect_classification_results()
        
        # Collect AMR detection results
        self._collect_amr_results()
        
        # Collect functional annotation results
        self._collect_annotation_results()
        
    def _collect_preprocessing_stats(self):
        """Collect preprocessing statistics"""
        qc_dir = self.input_dir / "01_preprocessing" / "qc"
        
        if qc_dir.exists():
            stats = []
            for stat_file in qc_dir.glob("*_nanostats.txt"):
                sample_stats = self._parse_nanostat_file(stat_file)
                if sample_stats:
                    stats.append(sample_stats)
            
            if stats:
                self.data['preprocessing']['stats'] = pd.DataFrame(stats)
                logger.info(f"Collected preprocessing stats for {len(stats)} samples")
    
    def _parse_nanostat_file(self, stat_file):
        """Parse NanoStat output file"""
        stats = {'sample': stat_file.stem.replace('_nanostats', '')}
        
        try:
            with open(stat_file, 'r') as f:
                for line in f:
                    if 'Mean read length:' in line:
                        stats['mean_length'] = float(line.split(':')[1].strip().replace(',', ''))
                    elif 'Median read length:' in line:
                        stats['median_length'] = float(line.split(':')[1].strip().replace(',', ''))
                    elif 'Number of reads:' in line:
                        stats['num_reads'] = int(line.split(':')[1].strip().replace(',', ''))
                    elif 'Read length N50:' in line:
                        stats['n50'] = float(line.split(':')[1].strip().replace(',', ''))
                    elif 'Total bases:' in line:
                        stats['total_bases'] = int(line.split(':')[1].strip().replace(',', ''))
                    elif 'Mean read quality:' in line:
                        stats['mean_quality'] = float(line.split(':')[1].strip())
        except Exception as e:
            logger.error(f"Error parsing {stat_file}: {e}")
            return None
        
        return stats
    
    def _collect_assembly_stats(self):
        """Collect assembly statistics"""
        assembly_dir = self.input_dir / "02_assembly" / "stats"
        
        if assembly_dir.exists():
            stats = []
            for stat_file in assembly_dir.glob("*_assemblystats.txt"):
                sample_stats = self._parse_assembly_stats(stat_file)
                if sample_stats:
                    stats.append(sample_stats)
            
            if stats:
                self.data['assembly']['stats'] = pd.DataFrame(stats)
                logger.info(f"Collected assembly stats for {len(stats)} samples")
    
    def _parse_assembly_stats(self, stat_file):
        """Parse assembly statistics file"""
        stats = {'sample': stat_file.stem.replace('_assemblystats', '')}
        
        try:
            with open(stat_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        metric, value = parts[0], parts[1]
                        if metric == 'sum':
                            stats['total_length'] = int(value)
                        elif metric == 'n':
                            stats['num_contigs'] = int(value)
                        elif metric == 'ave':
                            stats['avg_length'] = float(value)
                        elif metric == 'largest':
                            stats['largest_contig'] = int(value)
                        elif metric == 'N50':
                            stats['n50'] = int(value)
        except Exception as e:
            logger.error(f"Error parsing {stat_file}: {e}")
            return None
        
        return stats
    
    def _collect_classification_results(self):
        """Collect taxonomic classification results"""
        kraken_dir = self.input_dir / "03_classification" / "kraken2_reads"
        
        if kraken_dir.exists():
            taxa_abundance = defaultdict(lambda: defaultdict(float))
            
            for report_file in kraken_dir.glob("report_*.txt"):
                sample = report_file.stem.replace('report_', '')
                
                try:
                    with open(report_file, 'r') as f:
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 6:
                                percentage = float(parts[0])
                                rank = parts[3]
                                name = parts[5].strip()
                                
                                if rank in ['P', 'G', 'S']:  # Phylum, Genus, Species
                                    taxa_abundance[rank][(sample, name)] = percentage
                except Exception as e:
                    logger.error(f"Error parsing {report_file}: {e}")
            
            # Convert to DataFrames
            for rank, data in taxa_abundance.items():
                df_data = []
                for (sample, taxon), abundance in data.items():
                    df_data.append({
                        'sample': sample,
                        'taxon': taxon,
                        'abundance': abundance
                    })
                
                if df_data:
                    rank_name = {'P': 'phylum', 'G': 'genus', 'S': 'species'}[rank]
                    self.data['classification'][rank_name] = pd.DataFrame(df_data)
                    logger.info(f"Collected {rank_name} level classification data")
    
    def _collect_amr_results(self):
        """Collect AMR detection results"""
        amr_dir = self.input_dir / "05_amr_detection"
        
        # Collect ABRicate results
        abricate_dir = amr_dir / "abricate"
        if abricate_dir.exists():
            amr_genes = []
            
            for result_file in abricate_dir.glob("*_abricate.txt"):
                sample = result_file.stem.replace('_abricate', '')
                
                try:
                    df = pd.read_csv(result_file, sep='\t')
                    if not df.empty:
                        df['sample'] = sample
                        amr_genes.append(df)
                except Exception as e:
                    logger.error(f"Error reading {result_file}: {e}")
            
            if amr_genes:
                self.data['amr']['abricate'] = pd.concat(amr_genes, ignore_index=True)
                logger.info(f"Collected ABRicate results for {len(amr_genes)} samples")
    
    def _collect_annotation_results(self):
        """Collect functional annotation results"""
        eggnog_dir = self.input_dir / "04_annotation" / "eggnog"
        
        if eggnog_dir.exists():
            annotations = []
            
            for result_file in eggnog_dir.glob("*.emapper.annotations"):
                sample = result_file.stem.replace('.emapper', '')
                
                try:
                    df = pd.read_csv(result_file, sep='\t', comment='#')
                    if not df.empty:
                        df['sample'] = sample
                        annotations.append(df)
                except Exception as e:
                    logger.error(f"Error reading {result_file}: {e}")
            
            if annotations:
                self.data['annotation']['eggnog'] = pd.concat(annotations, ignore_index=True)
                logger.info(f"Collected eggNOG annotations for {len(annotations)} samples")
    
    def create_visualizations(self):
        """Create all visualizations for the report"""
        logger.info("Creating visualizations...")
        
        # Preprocessing visualizations
        if 'preprocessing' in self.data and 'stats' in self.data['preprocessing']:
            self._create_preprocessing_plots()
        
        # Assembly visualizations
        if 'assembly' in self.data and 'stats' in self.data['assembly']:
            self._create_assembly_plots()
        
        # Classification visualizations
        if 'classification' in self.data:
            self._create_classification_plots()
        
        # AMR visualizations
        if 'amr' in self.data:
            self._create_amr_plots()
    
    def _create_preprocessing_plots(self):
        """Create preprocessing visualization plots"""
        df = self.data['preprocessing']['stats']
        
        # Read statistics overview
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Read Count', 'Total Bases', 'Read Length Distribution', 'Quality Distribution')
        )
        
        # Read count bar plot
        fig.add_trace(
            go.Bar(x=df['sample'], y=df['num_reads'], name='Read Count'),
            row=1, col=1
        )
        
        # Total bases bar plot
        fig.add_trace(
            go.Bar(x=df['sample'], y=df['total_bases'], name='Total Bases'),
            row=1, col=2
        )
        
        # Read length box plot
        for sample in df['sample']:
            sample_data = df[df['sample'] == sample].iloc[0]
            fig.add_trace(
                go.Box(
                    y=[sample_data['mean_length']], 
                    name=sample,
                    showlegend=False
                ),
                row=2, col=1
            )
        
        # N50 bar plot
        fig.add_trace(
            go.Bar(x=df['sample'], y=df['n50'], name='N50'),
            row=2, col=2
        )
        
        fig.update_layout(height=800, showlegend=False)
        self.figures['preprocessing_stats'] = fig
    
    def _create_assembly_plots(self):
        """Create assembly visualization plots"""
        df = self.data['assembly']['stats']
        
        # Assembly statistics
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Number of Contigs', 'Total Assembly Length', 
                          'Contig N50', 'Largest Contig')
        )
        
        # Number of contigs
        fig.add_trace(
            go.Bar(x=df['sample'], y=df['num_contigs']),
            row=1, col=1
        )
        
        # Total length
        fig.add_trace(
            go.Bar(x=df['sample'], y=df['total_length']),
            row=1, col=2
        )
        
        # N50
        fig.add_trace(
            go.Bar(x=df['sample'], y=df['n50']),
            row=2, col=1
        )
        
        # Largest contig
        fig.add_trace(
            go.Bar(x=df['sample'], y=df['largest_contig']),
            row=2, col=2
        )
        
        fig.update_layout(height=800, showlegend=False)
        self.figures['assembly_stats'] = fig
    
    def _create_classification_plots(self):
        """Create taxonomic classification plots"""
        # Genus level abundance heatmap
        if 'genus' in self.data['classification']:
            df = self.data['classification']['genus']
            
            # Get top 20 most abundant genera
            top_genera = df.groupby('taxon')['abundance'].sum().nlargest(20).index
            df_top = df[df['taxon'].isin(top_genera)]
            
            # Pivot for heatmap
            heatmap_data = df_top.pivot(index='taxon', columns='sample', values='abundance').fillna(0)
            
            fig = go.Figure(data=go.Heatmap(
                z=heatmap_data.values,
                x=heatmap_data.columns,
                y=heatmap_data.index,
                colorscale='Viridis'
            ))
            
            fig.update_layout(
                title='Top 20 Genera Abundance Heatmap',
                xaxis_title='Sample',
                yaxis_title='Genus',
                height=600
            )
            
            self.figures['genus_heatmap'] = fig
    
    def _create_amr_plots(self):
        """Create AMR visualization plots"""
        if 'abricate' in self.data['amr']:
            df = self.data['amr']['abricate']
            
            # AMR genes per sample
            gene_counts = df.groupby('sample')['GENE'].count()
            
            fig = go.Figure(data=[
                go.Bar(x=gene_counts.index, y=gene_counts.values)
            ])
            
            fig.update_layout(
                title='AMR Genes Detected per Sample',
                xaxis_title='Sample',
                yaxis_title='Number of AMR Genes',
                height=400
            )
            
            self.figures['amr_counts'] = fig
    
    def generate_report(self):
        """Generate the final HTML report"""
        logger.info("Generating HTML report...")
        
        # HTML template
        template_str = '''
<!DOCTYPE html>
<html>
<head>
    <title>Nanopore Metagenomics Analysis Report</title>
    <meta charset="utf-8">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 40px;
            background-color: #f5f5f5;
        }
        .container {
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            max-width: 1200px;
            margin: 0 auto;
        }
        h1, h2, h3 {
            color: #333;
        }
        h1 {
            border-bottom: 3px solid #007bff;
            padding-bottom: 10px;
        }
        h2 {
            margin-top: 40px;
            border-bottom: 1px solid #ddd;
            padding-bottom: 5px;
        }
        .summary-box {
            background-color: #f8f9fa;
            border-left: 4px solid #007bff;
            padding: 15px;
            margin: 20px 0;
        }
        .metric {
            display: inline-block;
            margin: 10px 20px 10px 0;
        }
        .metric-value {
            font-size: 24px;
            font-weight: bold;
            color: #007bff;
        }
        .metric-label {
            color: #666;
            font-size: 14px;
        }
        table {
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }
        th {
            background-color: #007bff;
            color: white;
        }
        tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        .plot-container {
            margin: 20px 0;
        }
        .footer {
            margin-top: 50px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            text-align: center;
            color: #666;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Nanopore Metagenomics Analysis Report</h1>
        
        <div class="summary-box">
            <h3>Analysis Summary</h3>
            <div class="metric">
                <div class="metric-value">{{ n_samples }}</div>
                <div class="metric-label">Samples Analyzed</div>
            </div>
            <div class="metric">
                <div class="metric-value">{{ total_reads }}</div>
                <div class="metric-label">Total Reads</div>
            </div>
            <div class="metric">
                <div class="metric-value">{{ total_bases }}</div>
                <div class="metric-label">Total Bases</div>
            </div>
            <div class="metric">
                <div class="metric-value">{{ generation_date }}</div>
                <div class="metric-label">Report Generated</div>
            </div>
        </div>
        
        <h2>1. Preprocessing Results</h2>
        {% if preprocessing_table %}
        <h3>Read Statistics Summary</h3>
        {{ preprocessing_table }}
        {% endif %}
        
        {% if preprocessing_stats_plot %}
        <div class="plot-container">
            {{ preprocessing_stats_plot }}
        </div>
        {% endif %}
        
        <h2>2. Assembly Results</h2>
        {% if assembly_table %}
        <h3>Assembly Statistics Summary</h3>
        {{ assembly_table }}
        {% endif %}
        
        {% if assembly_stats_plot %}
        <div class="plot-container">
            {{ assembly_stats_plot }}
        </div>
        {% endif %}
        
        <h2>3. Taxonomic Classification</h2>
        {% if genus_heatmap_plot %}
        <h3>Genus-level Abundance</h3>
        <div class="plot-container">
            {{ genus_heatmap_plot }}
        </div>
        {% endif %}
        
        <h2>4. Antimicrobial Resistance Detection</h2>
        {% if amr_counts_plot %}
        <h3>AMR Genes per Sample</h3>
        <div class="plot-container">
            {{ amr_counts_plot }}
        </div>
        {% endif %}
        
        {% if amr_table %}
        <h3>AMR Genes Summary</h3>
        {{ amr_table }}
        {% endif %}
        
        <h2>5. Functional Annotation</h2>
        {% if annotation_summary %}
        {{ annotation_summary }}
        {% endif %}
        
        <div class="footer">
            <p>Generated by Nanopore Metagenomics Pipeline v1.0.0</p>
            <p>{{ generation_date }}</p>
        </div>
    </div>
</body>
</html>
        '''
        
        # Prepare template data
        template_data = {
            'generation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'n_samples': 0,
            'total_reads': 0,
            'total_bases': 0
        }
        
        # Add preprocessing data
        if 'preprocessing' in self.data and 'stats' in self.data['preprocessing']:
            df = self.data['preprocessing']['stats']
            template_data['n_samples'] = len(df)
            template_data['total_reads'] = f"{df['num_reads'].sum():,}"
            template_data['total_bases'] = f"{df['total_bases'].sum():,}"
            template_data['preprocessing_table'] = df.to_html(index=False, classes='table')
        
        # Add assembly data
        if 'assembly' in self.data and 'stats' in self.data['assembly']:
            df = self.data['assembly']['stats']
            template_data['assembly_table'] = df.to_html(index=False, classes='table')
        
        # Add AMR summary
        if 'amr' in self.data and 'abricate' in self.data['amr']:
            df = self.data['amr']['abricate']
            summary = df.groupby('sample')['GENE'].count().reset_index()
            summary.columns = ['Sample', 'AMR Genes']
            template_data['amr_table'] = summary.to_html(index=False, classes='table')
        
        # Add plots
        for plot_name, fig in self.figures.items():
            template_data[f'{plot_name}_plot'] = fig.to_html(include_plotlyjs=False, div_id=plot_name)
        
        # Render template
        template = Template(template_str)
        html_content = template.render(**template_data)
        
        # Write report
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(self.output_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Report generated: {self.output_file}")
    
    def run(self):
        """Run the complete report generation pipeline"""
        self.collect_data()
        self.create_visualizations()
        self.generate_report()


def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive metagenomics analysis report')
    parser.add_argument('--input-dir', required=True, help='Pipeline output directory')
    parser.add_argument('--output-file', required=True, help='Output HTML report file')
    parser.add_argument('--config', help='Configuration file (optional)')
    
    args = parser.parse_args()
    
    # Create report generator
    generator = MetagenomicsReportGenerator(args.input_dir, args.output_file)
    
    # Generate report
    generator.run()


if __name__ == '__main__':
    main()

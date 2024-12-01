# bactQC/cli.py

import click
from .core import Genome
import emoji
import os
from rich.console import Console
from rich.table import Table
from rich import box

# Initialize Rich console
console = Console()

# Define ASCII art as a raw multi-line string to prevent escape sequence warnings
ASCII_ART = r"""
   ___           __            _
  / _ )___ _____/ /____  ___  (_)__ _
 / _  / _ `/ __/ __/ _ \/ _ \/ / _ `/
/_____\_______/________/ .__/_/__,_/
/ __ \/ ___/ /_  __/___/___  / /__
/ /_/ / /__    / / / _ \/ _ \/ (_-<
\___\_\___/   /_/  \___/\___/_/___/

"""


@click.group()
def cli():
    """bactQC: A tool for bacterial genome quality control."""
    pass


@cli.command()
@click.option('--sample_name', help='Name of a sample to analyze')
@click.option('--input_dir', default='bactopia', type=click.Path(exists=True), help='Directory containing Bactopia outputs.')
@click.option('--min_primary_abundance', default=0.80, help='Minimum required abundance for the primary species.')
@click.option('--min_completeness', default=80, help='Minimum required completeness threshold.')
@click.option('--max_contamination', default=10, help='Maximum allowed contamination threshold.')
@click.option('--maximum_contigs', default=500, help='Maximum allowed number of contigs.')
@click.option('--minimum_n50', default=15000, help='Minimum required N50 contig length.')
@click.option('--min_q30_bases', default=0.90, help='Minimum required proportion of Q30 bases after filtering.')
@click.option('--min_coverage', default=30, help='Minimum required coverage after filtering.')
def run(sample_name, input_dir, min_primary_abundance, min_completeness, max_contamination, maximum_contigs, minimum_n50, min_q30_bases, min_coverage):
    """Run all quality control checks for a sample."""
    
    # Display ASCII Art
    console.print(ASCII_ART, style="bright_green bold")
    
    # Display Thank You and GitHub link
    console.print("Thanks for using Bactopia QC tools!", style="bright_green bold")
    console.print("Please report any issues on GitHub: https://github.com/maxlcummins/bactQC/issues", style="bright_green bold")
    
    # Initialize Genome object from core.py
    qc = Genome(sample_name, input_dir)
    try:
        console.print("Running QC analysis", style="bold green")
        qc.run(
            min_primary_abundance=min_primary_abundance,
            min_completeness=min_completeness,
            max_contamination=max_contamination,
            maximum_contigs=maximum_contigs,
            minimum_n50=minimum_n50,
            min_q30_bases=min_q30_bases,
            min_coverage=min_coverage
        )
        console.print("Determining if genome passed or failed", style="bold green")
        
        # Determine output filenames based on whether sample_name is provided
        if sample_name:
            output_prefix_results = f"{sample_name}_qc_results"
            output_prefix_thresholds = f"{sample_name}_qc_thresholds"
        else:
            output_prefix_results = "BactQC_results"
            output_prefix_thresholds = "BactQC_thresholds"
        
        results = qc.get_qc_results(output_prefix=output_prefix_results)
        thresholds = qc.get_qc_thresholds(output_prefix=output_prefix_thresholds)
    except Exception as e:
        console.print(f"Error during QC run: {e}", style="bold red")
        exit(1)
        
    # Get the absolute path of the input directory
    abs_input_path = os.path.abspath(input_dir)
    
    # Print target sample and input directory
    if sample_name:
        console.print(f"\n🦠🧬 Analysing Bactopia outputs for [bold]{sample_name}[/bold]")
    else:
        console.print("\n🦠🧬 Analysing Bactopia outputs for all samples")
    console.print(f"Checking [bold]{abs_input_path}[/bold]", style="bright_green")
    
    # Merge thresholds with results to include 'Detected species (Bracken)'
    if 'sample' in thresholds.columns and 'sample' in results.columns and 'Detected species (Bracken)' in results.columns:
        merged_thresholds = thresholds.merge(results[['sample', 'Detected species (Bracken)']], on='sample', how='left')
    else:
        console.print("Required columns for merging thresholds and results are missing.", style="bold red")
        exit(1)
    
    # Display QC Thresholds in Summarized Wide Format with Species-Specific Columns
    display_thresholds_summary(merged_thresholds)
    
    # Display QC Results
    display_qc_results(results)
    
    # Print the results file message correctly
    console.print(f"\n💾 Results written to [cyan]{output_prefix_results}.tsv[/cyan]", style="bold cyan")
    console.print(f"💾 Thresholds written to [cyan]{output_prefix_thresholds}.tsv[/cyan]", style="bold cyan")


def display_thresholds_summary(thresholds):
    """Display QC Thresholds in a summarized, well-formatted table with species-specific columns."""
    if thresholds.empty:
        console.print("No QC thresholds to display.", style="yellow")
        return
    
    console.print("\n❓ Quality Control Thresholds:", style="bold cyan")
    
    # Define columns to exclude
    columns_to_exclude = ['sample', 'mlst_expected_genus']
    # Determine columns to display
    display_columns = [col for col in thresholds.columns if col not in columns_to_exclude]
    
    if not display_columns:
        console.print("No QC thresholds to display after excluding specified columns.", style="yellow")
        return
    
    # Extract genus from 'Detected species (Bracken)'
    # Assuming 'Detected species (Bracken)' is in the format "Genus species", e.g., "Escherichia coli"
    thresholds['genus'] = thresholds['Detected species (Bracken)'].apply(
        lambda x: x.split()[0] if isinstance(x, str) and len(x.split()) >=1 else 'Unknown'
    )
    
    # Identify species-specific parameters
    # For simplicity, we assume 'assembly_scan_minimum_ungapped_length' and 'assembly_scan_maximum_ungapped_length' vary by species
    species_specific_params = [
        'assembly_scan_minimum_ungapped_length',
        'assembly_scan_maximum_ungapped_length'
    ]
    
    # Initialize a dictionary to hold parameters
    threshold_dict = {}
    
    # Handle non-species-specific parameters
    for col in display_columns:
        if col not in species_specific_params:
            # Assuming all non-species-specific parameters have the same value across samples
            unique_values = thresholds[col].unique()
            if len(unique_values) == 1:
                threshold_dict[col.replace('_', ' ').title()] = unique_values[0]
            else:
                # If values differ, it's unexpected; handle gracefully
                threshold_dict[col.replace('_', ' ').title()] = ', '.join(map(str, unique_values))
    
    # Handle species-specific parameters
    for param in species_specific_params:
        for sp in thresholds['genus'].unique():
            # Get the threshold value for this species and parameter
            species_threshold = thresholds[thresholds['genus'] == sp]
            if species_threshold.empty:
                continue  # No data for this species
            value = species_threshold[param].iloc[0]
            # Create a key with species name
            param_name = f"{param.replace('_', ' ').title()} ({sp})"
            threshold_dict[param_name] = value
    
    # Create a Rich table for thresholds summary
    thresholds_table = Table(title="", box=box.MINIMAL_DOUBLE_HEAD)
    thresholds_table.add_column("Parameter", style="bold green")
    thresholds_table.add_column("Value", style="cyan")
    
    for key, value in threshold_dict.items():
        thresholds_table.add_row(key, str(value))
    
    console.print(thresholds_table)


@cli.command()
@click.option('--sample_name', required=True, help='Name of a sample to analyze')
@click.option('--input_dir', default='bactopia', type=click.Path(exists=True), help='Directory containing Bactopia outputs.')
def get_assembly_size(sample_name, input_dir):
    """Retrieve total contig length from assembler results."""
    qc = Genome(sample_name, input_dir)
    try:
        qc.get_assembly_size(sample_name)
        assembly_size = qc.qc_data[sample_name].get('assembly_size', {})
    except Exception as e:
        console.print(f"Error retrieving assembly size: {e}", style="bold red")
        exit(1)
    
    if not assembly_size:
        console.print("No assembly size data available.", style="yellow")
        return
    
    # Create a Rich table for assembly size
    assembly_table = Table(title="Assembly Size", box=box.MINIMAL_DOUBLE_HEAD)
    assembly_table.add_column("Parameter", style="bold green")
    assembly_table.add_column("Value", style="cyan")
    
    for key, value in assembly_size.items():
        assembly_table.add_row(key.replace('_', ' ').title(), str(value))
    
    console.print(assembly_table)


@cli.command()
@click.option('--sample_name', required=True, help='Name of a sample to analyze')
@click.option('--input_dir', default='bactopia', type=click.Path(exists=True), help='Directory containing Bactopia outputs.')
@click.option('--min_primary_abundance', default=0.80, help='Minimum required abundance for the primary species.')
def check_bracken(sample_name, input_dir, min_primary_abundance):
    """Check Bracken results for a sample."""
    qc = Genome(sample_name, input_dir)
    try:
        qc.check_bracken(sample_name, min_primary_abundance)
        bracken_data = qc.qc_data[sample_name].get('bracken', {})
    except Exception as e:
        console.print(f"Error checking Bracken results: {e}", style="bold red")
        exit(1)
    
    if not bracken_data:
        console.print("No Bracken data available.", style="yellow")
        return
    
    # Create a Rich table for Bracken data
    bracken_table = Table(title="Bracken Data", box=box.MINIMAL_DOUBLE_HEAD)
    bracken_table.add_column("Parameter", style="bold green")
    bracken_table.add_column("Value", style="cyan")
    
    for key, value in bracken_data.items():
        # Format parameter name
        param_name = key.replace('_', ' ').title()
        bracken_table.add_row(param_name, str(value))
    
    console.print(bracken_table)


@cli.command()
@click.option('--sample_name', required=True, help='Name of a sample to analyze')
@click.option('--input_dir', default='bactopia', type=click.Path(exists=True), help='Directory containing Bactopia outputs.')
def check_mlst(sample_name, input_dir):
    """Check MLST results for a sample."""
    qc = Genome(sample_name, input_dir)
    try:
        qc.get_expected_genome_size(sample_name)
        expected_genus = qc.qc_data[sample_name]['genome_size']['organism_name'].split()[0]
        qc.check_mlst(sample_name, expected_genus)
        mlst_data = qc.qc_data[sample_name].get('mlst', {})
    except Exception as e:
        console.print(f"Error checking MLST results: {e}", style="bold red")
        exit(1)
    
    if not mlst_data:
        console.print("No MLST data available.", style="yellow")
        return
    
    # Create a Rich table for MLST data
    mlst_table = Table(title="MLST Data", box=box.MINIMAL_DOUBLE_HEAD)
    mlst_table.add_column("Parameter", style="bold green")
    mlst_table.add_column("Value", style="cyan")
    
    for key, value in mlst_data.items():
        # Format parameter name
        param_name = key.replace('_', ' ').title()
        mlst_table.add_row(param_name, str(value))
    
    console.print(mlst_table)


@cli.command()
@click.option('--sample_name', required=True, help='Name of a sample to analyze')
@click.option('--input_dir', default='bactopia', type=click.Path(exists=True), help='Directory containing Bactopia outputs.')
@click.option('--min_completeness', default=80, help='Minimum required completeness threshold.')
@click.option('--max_contamination', default=10, help='Maximum allowed contamination threshold.')
def check_checkm(sample_name, input_dir, min_completeness, max_contamination):
    """Check CheckM results for a sample."""
    qc = Genome(sample_name, input_dir)
    try:
        qc.check_checkm(sample_name, min_completeness, max_contamination)
        checkm_data = qc.qc_data[sample_name].get('checkm', {})
    except Exception as e:
        console.print(f"Error checking CheckM results: {e}", style="bold red")
        exit(1)
    
    if not checkm_data:
        console.print("No CheckM data available.", style="yellow")
        return
    
    # Create a Rich table for CheckM data
    checkm_table = Table(title="CheckM Data", box=box.MINIMAL_DOUBLE_HEAD)
    checkm_table.add_column("Parameter", style="bold green")
    checkm_table.add_column("Value", style="cyan")
    
    for key, value in checkm_data.items():
        # Format parameter name
        param_name = key.replace('_', ' ').title()
        checkm_table.add_row(param_name, str(value))
    
    console.print(checkm_table)


@cli.command()
@click.option('--sample_name', required=True, help='Name of a sample to analyze')
@click.option('--input_dir', default='bactopia', type=click.Path(exists=True), help='Directory containing Bactopia outputs.')
@click.option('--maximum_contigs', default=500, help='Maximum allowed number of contigs.')
@click.option('--minimum_n50', default=15000, help='Minimum required N50 contig length.')
def check_assembly_scan(sample_name, input_dir, maximum_contigs, minimum_n50):
    """Check assembly scan results for a sample."""
    qc = Genome(sample_name, input_dir)
    try:
        qc.check_assembly_scan(sample_name, maximum_contigs, minimum_n50)
        assembly_scan_data = qc.qc_data[sample_name].get('assembly_scan', {})
    except Exception as e:
        console.print(f"Error checking Assembly Scan results: {e}", style="bold red")
        exit(1)
    
    if not assembly_scan_data:
        console.print("No Assembly Scan data available.", style="yellow")
        return
    
    # Create a Rich table for Assembly Scan data
    assembly_scan_table = Table(title="Assembly Scan Data", box=box.MINIMAL_DOUBLE_HEAD)
    assembly_scan_table.add_column("Parameter", style="bold green")
    assembly_scan_table.add_column("Value", style="cyan")
    
    for key, value in assembly_scan_data.items():
        # Format parameter name
        param_name = key.replace('_', ' ').title()
        assembly_scan_table.add_row(param_name, str(value))
    
    console.print(assembly_scan_table)


@cli.command()
@click.option('--sample_name', required=True, help='Name of a sample to analyze')
@click.option('--input_dir', default='bactopia', type=click.Path(exists=True), help='Directory containing Bactopia outputs.')
@click.option('--min_q30_bases', default=0.90, help='Minimum required proportion of Q30 bases after filtering.')
@click.option('--min_coverage', default=30, help='Minimum required coverage after filtering.')
def check_fastp(sample_name, input_dir, min_q30_bases, min_coverage):
    """Check fastp quality control data for a sample."""
    qc = Genome(sample_name, input_dir)
    try:
        qc.check_fastp(sample_name, min_q30_bases, min_coverage)
        fastp_data = qc.qc_data[sample_name].get('fastp', {})
    except Exception as e:
        console.print(f"Error checking FastP data: {e}", style="bold red")
        exit(1)
    
    if not fastp_data:
        console.print("No FastP data available.", style="yellow")
        return
    
    # Create a Rich table for FastP data
    fastp_table = Table(title="FastP Data", box=box.MINIMAL_DOUBLE_HEAD)
    fastp_table.add_column("Parameter", style="bold green")
    fastp_table.add_column("Value", style="cyan")
    
    for key, value in fastp_data.items():
        # Format parameter name
        param_name = key.replace('_', ' ').title()
        fastp_table.add_row(param_name, str(value))
    
    console.print(fastp_table)


def display_thresholds_summary(thresholds):
    """Display QC Thresholds in a summarized, well-formatted table with species-specific columns."""
    if thresholds.empty:
        console.print("No QC thresholds to display.", style="yellow")
        return
    
    console.print("\n❓ Quality Control Thresholds:", style="bold cyan")
    
    # Define columns to exclude
    columns_to_exclude = ['sample', 'mlst_expected_genus']
    # Determine columns to display
    display_columns = [col for col in thresholds.columns if col not in columns_to_exclude]
    
    if not display_columns:
        console.print("No QC thresholds to display after excluding specified columns.", style="yellow")
        return
    
    # Extract genus from 'Detected species (Bracken)'
    # Assuming 'Detected species (Bracken)' is in the format "Genus species", e.g., "Escherichia coli"
    thresholds['genus'] = thresholds['Detected species (Bracken)'].apply(
        lambda x: x.split()[0] if isinstance(x, str) and len(x.split()) >=1 else 'Unknown'
    )
    
    # Identify species-specific parameters
    # For simplicity, we assume 'assembly_scan_minimum_ungapped_length' and 'assembly_scan_maximum_ungapped_length' vary by species
    species_specific_params = [
        'assembly_scan_minimum_ungapped_length',
        'assembly_scan_maximum_ungapped_length'
    ]
    
    # Initialize a dictionary to hold parameters
    threshold_dict = {}
    
    # Handle non-species-specific parameters
    for col in display_columns:
        if col not in species_specific_params:
            # Assuming all non-species-specific parameters have the same value across species
            unique_values = thresholds[col].unique()
            if len(unique_values) == 1:
                threshold_dict[col.replace('_', ' ').title()] = unique_values[0]
            else:
                # If values differ, list all unique values
                threshold_dict[col.replace('_', ' ').title()] = ', '.join(map(str, unique_values))
    
    # Handle species-specific parameters
    for param in species_specific_params:
        for sp in thresholds['genus'].unique():
            # Get the threshold value for this species and parameter
            species_threshold = thresholds[thresholds['genus'] == sp]
            if species_threshold.empty:
                continue  # No data for this species
            value = species_threshold[param].iloc[0]
            # Create a key with species name
            param_name = f"{param.replace('_', ' ').title()} ({sp})"
            threshold_dict[param_name] = value
    
    # Create a Rich table for thresholds summary
    thresholds_table = Table(title="", box=box.MINIMAL_DOUBLE_HEAD)
    thresholds_table.add_column("Parameter", style="bold green")
    thresholds_table.add_column("Value", style="cyan")
    
    for key, value in threshold_dict.items():
        thresholds_table.add_row(key, str(value))
    
    console.print(thresholds_table)


def display_qc_results(results):
    """Display QC Results in a well-formatted table with status and emojis."""
    if results.empty:
        console.print("No QC results to display.", style="yellow")
        return

    console.print("\n📊 Quality Control Results:", style="bold cyan")

    # Create a Rich table for results
    results_table = Table(title="", box=box.MINIMAL_DOUBLE_HEAD)

    # Add 'Sample' column separately
    results_table.add_column("Sample", style="bold magenta", no_wrap=True)

    # Exclude 'sample' from results.columns
    data_columns = [col for col in results.columns if col != 'sample']

    for col in data_columns:
        results_table.add_column(col.replace('_', ' ').title(), style="cyan", justify="center")

    # Iterate over each sample's results
    for idx, row in results.iterrows():
        sample_name = row.get('sample', 'N/A')
        row_values = [sample_name]  # Start with sample_name
        for col in data_columns:
            value = row[col]
            if isinstance(value, bool):
                if value:
                    status = emoji.emojize(":check_mark_button: [green]Passed[/green]")
                else:
                    status = emoji.emojize(":cross_mark: [red]Failed[/red]")
                row_values.append(status)
            else:
                row_values.append(str(value))
        results_table.add_row(*row_values)

    # Now, add a summary row at the end
    # Calculate the number of samples
    total_samples = len(results)

    # Initialize a summary dict
    summary = {'Sample': 'Total Passed'}

    for col in data_columns:
        # Only consider boolean columns
        if results[col].dtype == bool:
            num_passed = results[col].sum()
            percent_passed = (num_passed / total_samples) * 100
            summary[col.replace('_', ' ').title()] = f"{num_passed}/{total_samples} ({percent_passed:.1f}%)"
        else:
            summary[col.replace('_', ' ').title()] = ''  # Empty string for non-boolean columns

    # Build the summary row values
    summary_values = [summary['Sample']]
    for col in data_columns:
        summary_values.append(summary[col.replace('_', ' ').title()])

    # Add a separator row
    results_table.add_row(*['-' * 10] * len(summary_values))

    # Add the summary row
    results_table.add_row(*summary_values)

    console.print(results_table)


# Add commands to the CLI group
cli.add_command(run)
cli.add_command(get_assembly_size)
cli.add_command(check_bracken)
cli.add_command(check_mlst)
cli.add_command(check_checkm)
cli.add_command(check_assembly_scan)
cli.add_command(check_fastp)

def main():
    cli()

if __name__ == '__main__':
    main()

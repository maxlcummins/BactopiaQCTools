# bactQC/cli.py

import click
from .core import Genome
import emoji
import textwrap

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

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--taxid', default=None, help='Taxonomy ID of the species.')
def run(sample_name, input_dir, taxid):
    """Run all quality control checks for a sample."""
    qc = Genome(sample_name, input_dir, taxid)
    qc.run()
    results = qc.get_qc_results()
    thresholds = qc.get_qc_thresholds()
    
    # Show the ASCII art
    click.secho(ASCII_ART, fg='bright_green', bold=True)
    
    # Header
    click.secho("\n❓ Quality Control Thresholds:", fg='cyan', bold=True)
    for key, value in thresholds.items():
        click.echo(f"{key}: {value}")
        
    # Header
    click.secho("\n📊 Quality Control Results:", fg='cyan', bold=True)
    # Iterate over the results and display them with colors and emojis
    for check, passed in results.items():
        if check == 'sample':
            continue  # Skip the 'sample' key
        if passed:
            status = click.style('Passed', fg='green')
            emoji_icon = emoji.emojize(':check_mark_button:')
        else:
            status = click.style('Failed', fg='red')
            emoji_icon = emoji.emojize(':cross_mark:')
        click.echo(f"{check}: {status} {emoji_icon}")
        
    # Header
    click.secho(f"\n💾 Results written to {thresholds['sample']}_qc_results.tsv", fg='cyan', bold=True)

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--taxid', default=None, help='Taxonomy ID of the species.')
def get_expected_genome_size(sample_name, input_dir, taxid):
    """Retrieve expected genome size information."""
    qc = Genome(sample_name, input_dir, taxid)
    qc.get_expected_genome_size()
    click.echo("Expected genome size:")
    for key, value in qc.qc_data['genome_size'].items():
        click.echo(f"{key}: {value}")

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
def get_assembly_size(sample_name, input_dir):
    """Retrieve total contig length from assembler results."""
    qc = Genome(sample_name, input_dir)
    qc.get_assembly_size()
    click.echo("Assembly size:")
    for key, value in qc.qc_data['assembly_size'].items():
        click.echo(f"{key}: {value}")

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--min_primary_abundance', default=0.80, help='Minimum required abundance for the primary species.')
def check_bracken(sample_name, input_dir, min_primary_abundance):
    """Check Bracken results for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_bracken(min_primary_abundance)
    
    click.echo("Bracken data:")
    for key, value in qc.qc_data['bracken'].items():
        click.echo(f"{key}: {value}")

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.argument('expected_genus')
def check_mlst(sample_name, input_dir, expected_genus):
    """Check MLST results for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_mlst(expected_genus)
    click.echo("MLST data:")
    for key, value in qc.qc_data['mlst'].items():
        click.echo(f"{key}: {value}")

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--min_completeness', default=0.80, help='Minimum required completeness threshold.')
@click.option('--max_contamination', default=10, help='Maximum allowed contamination threshold.')
def check_checkm(sample_name, input_dir, min_completeness, max_contamination):
    """Check CheckM results for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_checkm(min_completeness, max_contamination)
    click.echo("CheckM data:")
    for key, value in qc.qc_data['checkm'].items():
        click.echo(f"{key}: {value}")

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--maximum_contigs', default=500, help='Maximum allowed number of contigs.')
@click.option('--minimum_n50', default=15000, help='Minimum required N50 contig length.')
def check_assembly_scan(sample_name, input_dir, maximum_contigs, minimum_n50):
    """Check assembly scan results for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_assembly_scan(maximum_contigs, minimum_n50)
    click.echo("Assembly-scan data:")
    for key, value in qc.qc_data['assembly_scan'].items():
        click.echo(f"{key}: {value}")

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--min_q30_bases', default=0.80, help='Minimum required proportion of Q30 bases after filtering.')
@click.option('--min_coverage', default=30, help='Minimum required coverage after filtering.')
def check_fastp(sample_name, input_dir, min_q30_bases, min_coverage):
    """Check fastp quality control data for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_fastp(min_q30_bases, min_coverage)
    click.echo("FastP data:")
    for key, value in qc.qc_data['fastp'].items():
        click.echo(f"{key}: {value}")

# Add commands to the CLI group
cli.add_command(run)
cli.add_command(get_expected_genome_size)
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

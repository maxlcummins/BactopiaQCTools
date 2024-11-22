# bactQC/cli.py

import click
from .core import Genome

@click.group()
def cli():
    """bactQC: A tool for bacterial genome quality control."""
    pass

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--taxid', default=None, help='Taxonomy ID of the species.')
def get_expected_genome_size(sample_name, input_dir, taxid):
    """Retrieve expected genome size information."""
    qc = Genome(sample_name, input_dir, taxid)
    qc.get_expected_genome_size()
    click.echo(qc.qc_results['genome_size'])

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
def get_assembly_size(sample_name, input_dir):
    """Retrieve total contig length from assembler results."""
    qc = Genome(sample_name, input_dir)
    qc.get_assembly_size()
    click.echo(qc.qc_results['assembly_size'])

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--min_primary_abundance', default=0.80, help='Minimum required abundance for the primary species.')
def check_bracken_results(sample_name, input_dir, min_primary_abundance):
    """Check Bracken results for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_bracken_results(min_primary_abundance)
    click.echo(qc.qc_results['bracken'])

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.argument('expected_genus')
def check_mlst_results(sample_name, input_dir, expected_genus):
    """Check MLST results for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_mlst_results(expected_genus)
    click.echo(qc.qc_results['mlst'])

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--min_completeness', default=0.80, help='Minimum required completeness threshold.')
@click.option('--max_contamination', default=10, help='Maximum allowed contamination threshold.')
def check_checkm_results(sample_name, input_dir, min_completeness, max_contamination):
    """Check CheckM results for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_checkm_results(min_completeness, max_contamination)
    click.echo(qc.qc_results['checkm'])

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--maximum_contigs', default=500, help='Maximum allowed number of contigs.')
@click.option('--minimum_n50', default=15000, help='Minimum required N50 contig length.')
def check_assembly_scan_results(sample_name, input_dir, maximum_contigs, minimum_n50):
    """Check assembly scan results for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_assembly_scan_results(maximum_contigs, minimum_n50)
    click.echo(qc.qc_results['assembly_scan'])

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--min_q30_bases', default=0.80, help='Minimum required proportion of Q30 bases after filtering.')
@click.option('--min_coverage', default=30, help='Minimum required coverage after filtering.')
def check_fastp_data(sample_name, input_dir, min_q30_bases, min_coverage):
    """Check fastp quality control data for a sample."""
    qc = Genome(sample_name, input_dir)
    qc.check_fastp_data(min_q30_bases, min_coverage)
    click.echo(qc.qc_results['fastp'])

@click.command()
@click.argument('sample_name')
@click.argument('input_dir')
@click.option('--taxid', default=None, help='Taxonomy ID of the species.')
def run(sample_name, input_dir, taxid):
    """Run all quality control checks for a sample."""
    qc = Genome(sample_name, input_dir, taxid)
    qc.run()
    results = qc.get_qc_results()
    for key, value in results.items():
        click.echo(f"{key}: {value}")

# Add commands to the CLI group
cli.add_command(get_expected_genome_size)
cli.add_command(get_assembly_size)
cli.add_command(check_bracken_results)
cli.add_command(check_mlst_results)
cli.add_command(check_checkm_results)
cli.add_command(check_assembly_scan_results)
cli.add_command(check_fastp_data)
cli.add_command(run)

def main():
    cli()

if __name__ == '__main__':
    main()

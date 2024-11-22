# bactQC/cli.py
import argparse
from bactQC import (
    get_expected_genome_size,
    get_assembly_size,
    check_bracken_results,
    check_mlst_results,
    check_checkm_results,
    check_assembly_scan_results,
    check_fastp_data,
    run,
)
import pprint

def main():
    parser = argparse.ArgumentParser(
        prog='bqc',
        description='bactQC Command-Line Interface',
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    # Subcommand: run
    parser_run = subparsers.add_parser('run', help='Run all QC checks')
    parser_run.add_argument('sample_name', type=str, help='Sample name')
    parser_run.add_argument('input_dir', type=str, help='Input directory path')
    parser_run.add_argument('--taxid', type=str, default=None, help='Taxonomy ID (optional)')

    # Subcommand: genome_size
    parser_genome_size = subparsers.add_parser('genome_size', help='Get expected genome size')
    parser_genome_size.add_argument('sample_name', type=str, help='Sample name')
    parser_genome_size.add_argument('input_dir', type=str, help='Input directory path')
    parser_genome_size.add_argument('--taxid', type=str, default=None, help='Taxonomy ID (optional)')

    # Subcommand: assembly_size
    parser_assembly_size = subparsers.add_parser('assembly_size', help='Get assembly size')
    parser_assembly_size.add_argument('sample_name', type=str, help='Sample name')
    parser_assembly_size.add_argument('input_dir', type=str, help='Input directory path')

    # Subcommand: check_bracken
    parser_bracken = subparsers.add_parser('check_bracken', help='Check Bracken results')
    parser_bracken.add_argument('sample_name', type=str, help='Sample name')
    parser_bracken.add_argument('input_dir', type=str, help='Input directory path')
    parser_bracken.add_argument('--min_primary_abundance', type=float, default=0.80, help='Minimum primary abundance (default: 0.80)')

    # Subcommand: check_mlst
    parser_mlst = subparsers.add_parser('check_mlst', help='Check MLST results')
    parser_mlst.add_argument('sample_name', type=str, help='Sample name')
    parser_mlst.add_argument('input_dir', type=str, help='Input directory path')
    parser_mlst.add_argument('expected_genus', type=str, help='Expected genus')

    # Subcommand: check_checkm
    parser_checkm = subparsers.add_parser('check_checkm', help='Check CheckM results')
    parser_checkm.add_argument('sample_name', type=str, help='Sample name')
    parser_checkm.add_argument('input_dir', type=str, help='Input directory path')
    parser_checkm.add_argument('--min_completeness', type=float, default=0.80, help='Minimum completeness (default: 0.80)')
    parser_checkm.add_argument('--max_contamination', type=float, default=10.0, help='Maximum contamination (default: 10.0)')

    # Subcommand: check_assembly_scan
    parser_assembly_scan = subparsers.add_parser('check_assembly_scan', help='Check assembly scan results')
    parser_assembly_scan.add_argument('sample_name', type=str, help='Sample name')
    parser_assembly_scan.add_argument('input_dir', type=str, help='Input directory path')
    parser_assembly_scan.add_argument('--maximum_contigs', type=int, default=500, help='Maximum contigs (default: 500)')
    parser_assembly_scan.add_argument('--minimum_N50', type=int, default=15000, help='Minimum N50 (default: 15000)')
    parser_assembly_scan.add_argument('--taxid', type=str, default=None, help='Taxonomy ID (optional)')

    # Subcommand: check_fastp
    parser_fastp = subparsers.add_parser('check_fastp', help='Check fastp data')
    parser_fastp.add_argument('sample_name', type=str, help='Sample name')
    parser_fastp.add_argument('input_dir', type=str, help='Input directory path')
    parser_fastp.add_argument('--min_q30_bases', type=float, default=0.80, help='Minimum Q30 bases (default: 0.80)')
    parser_fastp.add_argument('--min_coverage', type=int, default=30, help='Minimum coverage (default: 30)')

    args = parser.parse_args()
    pp = pprint.PrettyPrinter(indent=4)

    if args.command == 'run':
        result = run(args.sample_name, args.input_dir, args.taxid)
        print("QC Results:")
        pp.pprint(result)

    elif args.command == 'genome_size':
        result = get_expected_genome_size(args.sample_name, args.input_dir, args.taxid)
        print("Expected Genome Size Data:")
        pp.pprint(result)

    elif args.command == 'assembly_size':
        result = get_assembly_size(args.sample_name, args.input_dir)
        print("Assembly Size:")
        pp.pprint(result)

    elif args.command == 'check_bracken':
        result = check_bracken_results(args.sample_name, args.input_dir, args.min_primary_abundance)
        print("Bracken QC Results:")
        pp.pprint(result)

    elif args.command == 'check_mlst':
        result = check_mlst_results(args.sample_name, args.input_dir, args.expected_genus)
        print("MLST QC Results:")
        pp.pprint(result)

    elif args.command == 'check_checkm':
        result = check_checkm_results(
            args.sample_name,
            args.input_dir,
            min_completeness=args.min_completeness,
            max_contamination=args.max_contamination
        )
        print("CheckM QC Results:")
        pp.pprint(result)

    elif args.command == 'check_assembly_scan':
        result = check_assembly_scan_results(
            args.sample_name,
            args.input_dir,
            maximum_contigs=args.maximum_contigs,
            minimum_N50=args.minimum_N50,
            taxid=args.taxid
        )
        print("Assembly Scan QC Results:")
        pp.pprint(result)

    elif args.command == 'check_fastp':
        result = check_fastp_data(
            args.sample_name,
            args.input_dir,
            min_q30_bases=args.min_q30_bases,
            min_coverage=args.min_coverage
        )
        print("fastp QC Results:")
        pp.pprint(result)

    else:
        parser.print_help()

if __name__ == '__main__':
    main()

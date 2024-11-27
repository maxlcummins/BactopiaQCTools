# Bactopia QC Tools

Hello and welcome - note that all subcommands other than run are still under development. Don't bother trying to use them yet :)

## Installation

```bash
# Clone the repo
git clone https://github.com/maxlcummins/bactopiaQCtools.git

# Enter the directory
cd bactopiaQCtools

# Create an environment
mamba create -n bactQC python=3.12.6 -y

# Activate the environment
mamba activate bactQC

# Install bactopiaQCtools
pip install -e .
```

## Usage

After installing `bactopiaQCtools`, you can use it by running the following command:

```bash
bactQC --help
```

### Options

- `--input <path>`: Path to the input file or directory.
- `--output <path>`: Path to the output directory.
- `--threads <number>`: Number of threads to use.
- `--verbose`: Enable verbose output.

### Example

Here is an example of how to run `bactopiaQCtools`:

```bash
# For a Salmonella genome with the name Sample1
bactQC run Sample1 $(pwd)/bactopia_out --taxid 28901

# For a Salmonella genome with the name Sample1 with species autodetection
bactQC run Sample1 $(pwd)/bactopia_out
```

## Contributing

If you would like to contribute to `bactopiaQCtools`, please fork the repository and submit a pull request. We welcome all contributions!

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please open an issue on the [GitHub repository](https://github.com/maxlcummins/bactopiaQCtools).


# bactQC

![bactQC Logo](https://github.com/maxlcummins/bactQC/blob/main/assets/logo.png?raw=true)

**bactQC** is a command-line tool designed for downstream processing of analysis performed by [Bactopia](https://github.com/bactopia/bactopia).

It integrates multiple QC checks, including Bracken, MLST, CheckM, Assembly Scan, and Fastp, to ensure high-quality genomic data for downstream analyses.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Run All QC Checks](#run-all-qc-checks)
  - [Individual QC Commands](#individual-qc-commands)
    - [Get Expected Genome Size](#get-expected-genome-size)
    - [Get Assembly Size](#get-assembly-size)
    - [Check Bracken](#check-bracken)
    - [Check MLST](#check-mlst)
    - [Check CheckM](#check-checkm)
    - [Check Assembly Scan](#check-assembly-scan)
    - [Check Fastp](#check-fastp)
- [Command Reference](#command-reference)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)



## Features

- **Comprehensive QC Checks**: Integrates Bracken, MLST, CheckM, Assembly Scan, and Fastp for thorough quality assessment.
- **Flexible Genus Specification**: Determines the sample genus from Bactopia mash outputs or allows manual specification of an expected genus.
- **Rich Output Formatting**: Utilizes the Rich library for enhanced, visually appealing console outputs.
- **Error Handling**: Provides clear and informative error messages to guide users.
- **Modular Design**: Easily extendable to incorporate additional QC checks in the future.



## Installation

### Prerequisites

- **Python 3.7 or higher**: Ensure Python is installed. You can download it from the [official website](https://www.python.org/downloads/).
- **pip**: Python package installer. It usually comes with Python.

### Clone the Repository

```bash
# Clone our repo
git clone https://github.com/maxlcummins/bactQC.git

# Enter it
cd bactQC
```

### Install Dependencies
It's recommended to use a virtual environment to manage dependencies.

```bash
# Create a environment with mamba
mamba create -n bactQC python>=3.6

# Activate the environment
mamba activate bactQC

# Install required packages - Ensure you're in the bactQC directory containing setup.py
pip install .
```

Alternatively, you can install the dependencies manually

```bash
pip install click rich emoji pandas requests
```
###  Usage

There are two approaches to running this tools.

Firstly, the user can analyse a given sample based on any specific tool by executing a `bactQC` subcommand:

```
bactQC [OPTIONS] COMMAND [ARGS]
```

bactQC commands include `check-assembly-scan`, `check-bracken`, `check-checkm`, `check-fastp`, `check-mlst`, `get-assembly-size`, and `get-expected-genome-size`.

This will return key columns from the tool of interest as well as whether the genome passed QC for that tool (based on either default QC thresholds or user defined thresholds with the available flags).

For example, running `check-mlst` will return the following:

```
$bactQC check-mlst Genome123 /home/username/path/to/bactopia_output_directory

INFO:bactQC.core:Derived expected_genus from genome_size after fetching: Salmonella
INFO:bactQC.core:Loaded cached scheme_species_map.tab
INFO:bactQC.core:MLST scheme matches the expected genus schemes.
              MLST Data              
              ╷                      
  Parameter   │ Value                
 ═════════════╪═════════════════════ 
  Scheme      │ senterica_achtman_2  
  St          │ 3774                 
  Allele1     │ aroC(147)            
  Allele2     │ dnaN(14)             
  Allele3     │ hemD(21)             
  Allele4     │ hisD(88)             
  Allele5     │ purE(6)              
  Allele6     │ sucA(676)            
  Allele7     │ thrA(12)             
  Passed mlst │ True 
```

Secondly is via the `run` subcommand.This is intended to tell the user if the genome has passed or failed each of the quality control analyses.

This will be printed to the screen along with the quality control thresholds utilised:

```
$bactQC run Genome123 bactopia_output_directory

🦠🧬 Analysing Bactopia outputs for Genome123
Checking /home/username/path/to/bactopia_output_directory

❓ Quality Control Thresholds:
                ╷                                   
  QC Check      │ Parameters                        
 ═══════════════╪══════════════════════════════════ 
  Bracken       │ Min primary abundance: 0.8        
  Mlst          │ Expected genus: Salmonella        
  Checkm        │ Max contamination: 10             
                │ Min completeness: 80              
  Assembly scan │ Maximum contigs: 500              
                │ Minimum n50: 15000                
                │ Minimum ungapped length: 4100000  
                │ Maximum ungapped length: 6000000  
  Fastp         │ Min q30 bases: 0.9                
                │ Min coverage: 30                  
                ╵                                   

📊 Quality Control Results:
                ╷            
  QC Check      │  Status    
 ═══════════════╪═══════════ 
  Bracken       │ Passed ✅  
  Mlst          │ Passed ✅  
  Checkm        │ Passed ✅  
  Assembly scan │ Passed ✅  
  Fastp         │ Passed ✅  
                ╵            

💾 Results written to Genome123_qc_results.tsv
```

### Outputs

By default, it will also write to file a tab-delimited file containing simlar to the following:

|sample   |Detected species (Bracken)|Detected species (Mash)|bracken|mlst|checkm|assembly_scan|fastp|
|:--------|:-------------------------|:----------------------|:------|:---|:-----|:------------|:----|
|Genome123|Salmonella enterica       |Salmonella enterica    |True   |True|True  |True         | True|

The columns should be interpreted as follows:

| Column                   | Meaning                                               |
|:-------------------------|:------------------------------------------------------|
|Detected species (Bracken)|Which species was detected by Bactopia's run of Bracken|
|Detected species (Mash)   |Which species was detected by Bactopia's run of Mash   |
|bracken                   |Whether bracken QC passed                              |
|mlst                      |Whether mlst QC passed                                 |
|checkm                    |Whether checkm QC passed                               |
|assembly_scan             |Whether assembly_scan QC passed                        |
|fastp                     |Whether fastp QC passed                                |

### Examples

1. Running All QC Checks
```bash
bactQC run Genome_1 /data/genomes/Salmonella --taxid 28901
```

2. Performing MLST Check with Derived Genus
```bash
bactQC check-mlst Genome_1 /data/genomes/Salmonella
```

3. Performing MLST Check with Specified Genus
```bash
bactQC check-mlst Genome_1 /data/genomes/Salmonella --expected_genus Salmonella
```

4. Checking Fastp Quality, specifying minimum q30 score and coverage (defaults are 0.9 and 30)
```bash
bactQC check-fastp Genome_1 /data/genomes/Salmonella --min_q30_bases 0.95 --min_coverage 40
```

### Contributing
Contributions are welcome! To contribute to bactQC, please follow these guidelines:

### License
This project is licensed under the MIT License.

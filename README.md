# SV Bed Converter

This Python script converts a tab-separated input file containing structural variation (SV) data into a BED format file with additional information. It categorizes SVs based on their type and generates separate BED files for non-translocations and inter-chromosomal translocations.

## Usage 

To use this script, you need to provide two command-line arguments:
1. **Input file path**: The path to the input SV data file.
2. **Output file path**: The path where you want to save the generated BED file.

Here's an example of how to run the script:

```bash
python script_name.py input_filepath output_filepath
```

Replace script_name.py with the actual name of your Python script, input_filepath with the path to your input file, and output_filepath with the desired output file path.

## Input File
The input file should be a tab-separated text file containing the following columns:

* chrom1: Chromosome 1.
* position1: Start position on chromosome 1.
* position2: End position on chromosome 1.
* size: Size of the SV.
* type: Type of SV.
* classification: Classification of the SV.
* replicate: Replicate information.
Other columns (optional): Additional columns with data related to the SV.

## Output File
The script generates a BED format output file with the following columns:

* chrom: Chromosome.
* chromStart: Start position on the chromosome.
* chromEnd: End position on the chromosome.
* name: Custom name for the SV.
* score: Score (empty).
* strand: Strand (empty).
* blockCount: Block count (empty).
* blockSizes: Block sizes (empty).
* color: Color based on SV type.
The output file will be saved in the specified output file path.

## Dependencies
This script requires the following Python packages:

```
sys
pandas
os
pathlib
```
Make sure you have these packages installed in your Python environment before running the script.

## Support
This tool is maintained by Syukri Shukor. Please email Syukri Shukor, Alex Chitsazan, or Andy Pang for any questions.
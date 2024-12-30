# Transcript to Genomic Coordinates Conversion

The objective is to write software that translates transcript coordinates to genomic coordinates. To achieve this goal efficiently, we implement a data structure that stores the transcript coordinates and the corresponding genomic coordinates. Here we use [interval trees](https://en.wikipedia.org/wiki/Interval_tree) to efficiently query the genomic coordinates for a given transcript coordinate. The interval tree data structure allows us to quickly find the genomic coordinates that overlap with a given transcript coordinate. We implement a function that translates transcript coordinates to genomic coordinates using interval trees. The function takes a list of transcript coordinates as input and returns the corresponding genomic coordinates. We test the function with a list of transcript coordinates and verify that the output matches the expected genomic coordinates. The implementation is efficient and can handle large datasets with ease.

## Key Features
- **Coordinate Conversion**: Translate transcript coordinates to genomic coordinates.
- **Efficient Implementation**: Utilizes an interval tree data structure for fast lookups.
- **Scalability**: Handles large datasets with ease.


## Input File Format

### Transcript File
A tab-separated file with no header. Columns:
1. **transcript name** (string): Identifier for the transcript.
2. **chromosome name** (string): Associated chromosome.
3. **genomic start** (integer): Starting position on the genome.
4. **CIGAR string** (string): Encodes alignment details.

### Query File
A tab-separated file with no header. Columns:
1. **transcript name** (string): Identifier for the transcript.
2. **transcript position** (integer): Position on the transcript.

---

## Output File Format
A tab-separated file with no header. Columns:
1. **transcript name** (string): Identifier for the transcript.
2. **transcript position** (integer): Position on the transcript.
3. **chromosome name** (string): Associated chromosome.
4. **genomic position** (integer): Position on the genome.

## Requirements
- Python 3.6 or higher
- intervaltree
- pandas

## Usage

To run the script, use the following command:

```bash
cd src/
python main.py -t ../data/transcripts.tsv -q ../data/queries.tsv -o ../output/output.tsv
```
For a stdout visualization of the CIGAR string, see the following example:
```bash
python main.py -t ../data/transcripts.tsv -q ../data/queries.tsv -o ../output/output.tsv -d
```


## Testing
To test the code, I created unit tests to test the transcript class. Special attention was given to parsing the CIGAR and ensuring the correct transcript and genomic start values. The tests were run using unittest module. The tests cover the following scenarios:
- Test the transcript class with a simple example
- Test the transcript class with cases that could be problematic
- Test to see genomic values were positive integers
- Test to see if the transcript values were positive integers and within bounds of the mapped transcript

In addition, tests were run to ensure the user could catch errors in input files with missing values or incorrect values. These input files can be found in the test_input folder. One could go a step further and integrate test to check the input file format into the unittest test suite.

To run the tests, execute the following command:
```bash
cd src/
python -m unittest test_transcript.py
```

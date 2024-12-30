# Transcript to Genomic Coordinates Conversion

The objective is to write software that translates transcript coordinates to genomic coordinates. To achieve this goal efficiently, we implement a data structure that stores the transcript coordinates and the corresponding genomic coordinates. Here we use [interval trees](https://en.wikipedia.org/wiki/Interval_tree) to efficiently query the genomic coordinates for a given transcript coordinate. The simplest form of an interval tree can be though of as a binary search tree. This data structure allows you to query a value in binary search tree efficiently.  That is because binary tree search lookup structures the nodes in such a way that each comparison skips half the reminaing tree. In its simplest form, the time complexity of a binary search tree query is *O(log(n))* for *n* nodes. Interval trees are a more complex form of binary search trees that allow you to store intervals and query for overlapping intervals. Because of the added complexity, query time is *O(log(n) + m)* where *n* is the number of intervals in the collection and *m* is the number of intervals produce by the query. Creation of the interval tree data structure is *O(nlog(n))* This is useful for our problem because we can store the genomic coordinates as intervals and query for the genomic coordinates that overlap with a given transcript coordinate. 

Creation of the interval tree is a one time cost usually so this data structure would be better suited to larger datasets. In this implementation, we create interval trees for each transcript mapping which is why I chose to create a custom transcript class.  A more efficient solution utlizing the strengths of binary/iterval tree search would be to group transcripts together along their genomic coordinates during creation of the interval tree. Nonetheless, one easy improvement using the current implementation where each transcript considered separately would be to group queries by transcript.  This would might reduce the time it takes to retrieve the genomic coordinates for each transcript.

Here we implement functions that translates transcript coordinates to genomic coordinates using interval trees. The functions take transcript coordinates as input and returns the corresponding genomic coordinates.. The implementation is efficient and can handle large datasets with ease.

## Key Features
- **Coordinate Conversion**: Translate transcript coordinates to genomic coordinates.
- **Efficient Implementation**: Utilizes an interval tree data structure for fast lookups.
- **Scalability**: Handles large datasets with ease.

## Assumptions
- Transcripts are always mapped from genomic 5’ to 3’.
- The CIGAR string is in the format of the SAM file format.
- I am not really dealing with the S, P, H options in the CIGAR string. 


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
python main.py --transcript_file ../data/transcripts.tsv --query_file ../data/queries.tsv --output_file ../output/output.tsv
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

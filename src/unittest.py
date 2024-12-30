import unittest
import sys
sys.path.append("./src")
from transcript import Transcript

def test_cases():
    # Test cases
    # TODO create more test cases
    list_of_transcripts = [
        Transcript(
            transcript_id="TR1",
            chromosome="CHR1",
            genomic_start=3,
            cigar="8M7D6M2I2M11D7M",
        ),
        Transcript(
            transcript_id="TR2", chromosome="CHR2", genomic_start=10, cigar="20M"
        )
    ]
    return list_of_transcripts

def test_valid_cigar(transcript):
    # Check that the cigar
    valid_cigar = transcript.valid_cigar()
    assert valid_cigar == True
    transcript.visualize_cigar()

def test_coordinate(transcript, tx_start, expectd_genomic_coordinate):
    # Check that the genomic coordinate is correct
    transcript.precomputeIntervals2()
    coord = transcript.get_genomic_coordinates_from_precomputed_intervals2(tx_start)
    print(f"Genomic coordinate: {coord}")
    assert coord == expectd_genomic_coordinate

def main():
    # TODO maybe restructure so I don't have to loop through transcripts twice. 
    transcript_start = [4, 0]
    genomic_start = [7, 10]
    for i, transcript in enumerate(test_cases()):
        test_valid_cigar(transcript)
        test_coordinate(transcript, transcript_start[i], genomic_start[i])
    transcript_start = [13, 10]
    genomic_start = [23, 20]
    for i, transcript in enumerate(test_cases()):
        test_valid_cigar(transcript)
        test_coordinate(transcript, transcript_start[i], genomic_start[i])

if __name__ == "__main__":
    unittest.main()

import unittest
import sys
sys.path.append("..")  # Adjust the path to point to the parent directory

from transcript import Transcript

class TestTranscript(unittest.TestCase):
    def setUp(self):
        self.valid_transcripts = [
            Transcript(
                transcript_id="TR1",
                chromosome="CHR1",
                genomic_start=3,
                cigar="8M7D6M2I2M11D7M",
            ),
            Transcript(
                transcript_id="TR2", 
                chromosome="CHR2", 
                genomic_start=10, 
                cigar="20M"
            )
        ]
        self.invalid_transcripts = [
            {"transcript_id": "TR3", "chromosome": "CHR1", "genomic_start": "3", "cigar": "8Y7D6M2I2M11D7M"},
            {"transcript_id": "TR4", "chromosome": "CHR2", "genomic_start": -10, "cigar": "M"}
        ]

    def test_valid_cigar(self):
        for transcript in self.valid_transcripts:
            self.assertTrue(transcript.parse_cigar())
            transcript.visualize_cigar()

    def test_invalid_cigar(self):
        for invalid_transcript in self.invalid_transcripts:
            with self.assertRaises(ValueError):
                transcript = Transcript(**invalid_transcript)
                self.assertFalse(transcript.parse_cigar())


    def test_coordinate(self):
        transcript_start = [4, 0]
        genomic_start = [7, 10]
        for i, transcript in enumerate(self.valid_transcripts):
            transcript.precomputeIntervals2()
            coord = transcript.get_genomic_coordinates_from_precomputed_intervals2(transcript_start[i])
            self.assertEqual(coord, genomic_start[i])
    
    def test_coordinate_negative_tx_start(self):
        for i, transcript in enumerate(self.valid_transcripts):
            with self.assertRaises(ValueError):
                transcript.get_genomic_coordinates_from_precomputed_intervals2(-1)

    def test_bad_cases(self):
        for invalid_transcript in self.invalid_transcripts:
            with self.assertRaises(ValueError):
                Transcript(**invalid_transcript)


if __name__ == "__main__":
    unittest.main()
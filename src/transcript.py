"""
This class represents a transcript. We initialize it with basic information about the transcript, such as its ID, 
chromosome, genomic start position, and CIGAR string. 
"""


class Transcript:
    def __init__(self, transcript_id, chromosome, genomic_start, cigar):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.genomic_start = genomic_start
        self.cigar = cigar  # CIGAR string

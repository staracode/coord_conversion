"""
This class represents a transcript. We initialize it with basic information about the transcript, such as its ID, 
chromosome, genomic start position, and CIGAR string. 
"""

import re
class Transcript:
    def __init__(self, transcript_id, chromosome, genomic_start, cigar):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.genomic_start = int(genomic_start)
        self.cigar = cigar  # CIGAR string
        self.strand = "+"  # Assume all transcripts are on the positive strand
        
    
    def get_genomic_coordinate(self, tx_start):
        """
        Given a transcript coordinate (tx_start), return the corresponding genomic coordinate.
        """
        genomic_pos = self.genomic_start
        for length, operation in self.parse_cigar():
            length = int(length)
            if tx_start <= 0:
                return genomic_pos + tx_start
            if operation in ["M", "="]:
                genomic_pos += length
                tx_start -= length
            elif operation in ["D", "N"]:
                genomic_pos += length
            elif operation in ["I", "S"]: #what is S?
                tx_start -= length
            else:   # Skipped regions
                pass
        return genomic_pos  + tx_start

    def parse_cigar(self):
        """
        Parse the CIGAR string and return a list of tuples (length, operation).
        """
        cigar_tuples = re.findall(r'([0-9]+)([MIDNX=]{1})', self.cigar)
        return cigar_tuples

"""
This class represents a transcript. We initialize it with basic information about the transcript, such as its ID, 
chromosome, genomic start position, and CIGAR string. 
"""

import re
from intervaltree import Interval, IntervalTree

class Transcript:
    def __init__(self, transcript_id, chromosome, genomic_start, cigar):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.genomic_start = int(genomic_start)
        self.cigar = cigar
        self.cigar_tuples = self.parse_cigar()  # CIGAR string
        self.strand = "+"  # Assume all transcripts are on the positive strand
        self.tx_intervals = []
        self.gx_intervals = []
        self.exonIntervals = None

    def valid_cigar(self):
        """
        Check if the CIGAR string is valid.
        """
        cigar_operations = set(["M", "I", "D", "N", "S", "=", "X"])
        for length, operation in self.parse_cigar():
            if operation not in cigar_operations:
                return False
            if length <= 0:
                return False
        return True

    def parse_cigar(self):
        """
        Parse the CIGAR string and return a list of tuples (length, operation).
        """
        cigar_tuples = re.findall(r"([0-9]+)([MIDNX=]{1})", self.cigar)
        cigar_tuples = [(int(length), operation) for length, operation in cigar_tuples]
        return cigar_tuples

    def get_genomic_coordinate(self, tx_start):
        """
        Given a transcript coordinate (tx_start), return the corresponding genomic coordinate.
        """
        genomic_pos = self.genomic_start
        for length, operation in self.cigar_tuples:
            if tx_start <= 0:
                return genomic_pos + tx_start
            if operation in ["M", "=", "X"]:
                genomic_pos += length
                tx_start -= length
            elif operation in ["D", "N"]:
                genomic_pos += length
            elif operation in ["I", "S"]:  # what is S?
                tx_start -= length
            else:  # Skipped regions
                pass
        # TODO: if tx_start exceeds the length of the transcript, return the genomic position?
        return genomic_pos + tx_start

    def visualize_cigar(self):
        """
        Visualize the CIGAR string.
        """
        (genomic_string, query_string, count_string) = ("", "", "")
        for length, operation in self.cigar_tuples:
            spaces = (" " * (length - len(str(length)))) + str(
                length
            )  # Add spaces for alignment; won't align if length is less than the size of str(length)
            count_string += spaces
            if operation in ["M", "=", "X"]:
                genomic_string += "M" * int(length)
                query_string += "M" * int(length)
            elif operation in ["D", "N"]:
                genomic_string += "D" * int(length)
                query_string += "-" * int(length)
            elif operation in ["I", "S"]:
                genomic_string += "-" * int(length)
                query_string += "I" * int(length)
            else:
                pass
        print(f"Genomic: {genomic_string}")
        print(f"Query:   {query_string}")
        print(f"Count:   {count_string}")

    def precompute_tx_intervals(self):
        """
        Precompute the intervals of the transcript.
        """
        tx_intervals = []
        gx_intervals = []
        tx_pos = 0
        gx_pos = self.genomic_start
        for length, operation in self.cigar_tuples:
            if operation in ["M", "=", "X"]:
                tx_intervals.append((tx_pos, tx_pos + length - 1))
                gx_intervals.append((gx_pos, gx_pos + length - 1))
                tx_pos += length
                gx_pos += length
            elif operation in ["D", "N"]:
                tx_intervals.append(
                    (-tx_pos, -tx_pos)
                )  # negative values to indicate skipped regions
                gx_intervals.append((gx_pos, gx_pos + length - 1))
                gx_pos += length
            elif operation in ["I", "S"]:
                tx_intervals.append((tx_pos, tx_pos + length - 1))
                gx_intervals.append(
                    (-gx_pos, -gx_pos)
                )  # negative values to indicate skipped regions
                tx_pos += length
            else:
                pass

        print(f"Genomic intervals: {gx_intervals}")
        print(f"Transcript intervals: {tx_intervals}")
        # TODO: explore data structures to store intervals
        self.tx_intervals = tx_intervals
        self.gx_intervals = gx_intervals

    def get_genomic_coordinates_from_precomputed_intervals(self, tx_start):
        """
        Given a transcript coordinate (tx_start), return the corresponding genomic coordinate using precomputed intervals.
        """
        for i, tx_interval in enumerate(self.tx_intervals):
            tx_interval_start = tx_interval[0]
            tx_interval_end = tx_interval[1]
            if tx_start >= tx_interval_start and tx_start <= tx_interval_end:
                tx_pos = tx_start - tx_interval_start
                genomic_pos = self.gx_intervals[i][0] + tx_pos
                return genomic_pos

    def precomputeIntervals2(self):
        """
        Precompute the intervals of the transcript using the intervaltree library.
        """
        exonIntervals = IntervalTree()
        tx_pos = 0
        gx_pos = self.genomic_start
        for length, operation in self.cigar_tuples:
            if operation in ["M", "=", "X"]:
                exonIntervals.addi(tx_pos, tx_pos + length, gx_pos)
                tx_pos += length
                gx_pos += length
            elif operation in ["D", "N"]:  # what is N?
                #exonIntervals.addi(tx_pos, tx_pos, gx_pos + length)
                gx_pos += length
            elif operation in ["I", "S"]:
                tx_pos += length
            else:       
                pass
        self.exonIntervals = exonIntervals
    
    def get_genomic_coordinates_from_precomputed_intervals2(self, tx_start):
        """
        Given a transcript coordinate (tx_start), return the corresponding genomic coordinate using precomputed intervals.
        """
        for interval in self.exonIntervals[tx_start]:
            return interval.data + (tx_start - interval.begin)
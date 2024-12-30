"""
This class represents a transcript. We initialize it with basic information about the transcript, such as its ID, 
chromosome, genomic start position, and CIGAR string. 
"""

import re
from intervaltree import Interval, IntervalTree


class Transcript:
    def __init__(
        self, transcript_id: str, chromosome: str, genomic_start: int, cigar: str
    ):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        if not isinstance(genomic_start, int):
            raise ValueError(
                "genomic_start {} must be an integer".format(genomic_start)
            )
        if genomic_start < 0:
            raise ValueError(
                "genomic_start {} must be a non-negative integer".format(genomic_start)
            )
        self.genomic_start = genomic_start
        self.cigar = cigar.lstrip().rstrip() # Remove trailing whitespace
        self.cigar_tuples = self.parse_cigar()  # List of tuples (length, operation)
        self.strand = (
            "+"  # Assume that the transcript is always mapped from genomic 5’ to 3’.
        )
        # self.tx_intervals = []  # List of tuples (start, end)
        # self.gx_intervals = []  # List of tuples (start, end)
        self.exonIntervals = None
        self.tx_len = 0

    def valid_cigar(self, cigar_tuple):
        """
        Check if the CIGAR string is valid.
        """
        length, operation = cigar_tuple
        cigar_operations = set(["M", "I", "D", "N", "S", "H", "P", "=", "X"])
        if operation not in cigar_operations:
            raise ValueError(f"Invalid operation {operation} in CIGAR string: double CIGAR cigar {self.cigar}")
        if (
            length <= 0
        ):  # I think technically zero is allowed but I can't understand the logic behind it
            raise ValueError(f"Invalid length {length} in CIGAR string: double check CIGAR {self.cigar}")

    def parse_cigar(self):
        """
        Parse the CIGAR string and return a list of tuples (length, operation).
        """
        cigar_tuples = re.findall(
            r"([0-9]+)([A-Z|a-z|:/?#@!$&" "'" "()*+,;=%[]{1})", self.cigar
        )
        cigar_tuples_check = [
            (int(length), operation) for length, operation in cigar_tuples
        ]
        # check if the cigar string is valid or raise an error
        if cigar_tuples_check == []:
            raise ValueError("Invalid CIGAR string {}".format(self.cigar))
        new_cigar = ""
        for one_cigar_tuple in cigar_tuples_check:
            self.valid_cigar(one_cigar_tuple)
            new_cigar += str(one_cigar_tuple[0]) + str(one_cigar_tuple[1])
        # check if the new cigar string is the same as the original cigar string
        # if regex is not finding the cigar tuple pair, it drops the cigar tuple pair
        # this will cause the new_cigar to be different from the original cigar
        #print (new_cigar)
        if len(new_cigar) != len(self.cigar):
            raise ValueError(
                f"HereInvalid CIGAR string {self.cigar}: double check CIGAR {new_cigar}"
            )
        return cigar_tuples_check

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
                pass  # I don't want to deal with  H, P because they don't affect the genomic position
        print(f"Genomic: {genomic_string}")
        print(f"Query:   {query_string}")
        print(f"Count:   {count_string}")

    # def precompute_tx_intervals(self):
    #     """
    #     Precompute the intervals of the transcript.
    #     I don't use this because I don't consider it as effcient but it helps me check my logic.
    #     """
    #     tx_intervals = []
    #     gx_intervals = []
    #     tx_pos = 0
    #     gx_pos = self.genomic_start
    #     for length, operation in self.cigar_tuples:
    #         if operation in ["M", "=", "X"]:
    #             tx_intervals.append((tx_pos, tx_pos + length - 1))
    #             gx_intervals.append((gx_pos, gx_pos + length - 1))
    #             tx_pos += length
    #             gx_pos += length
    #         elif operation in ["D", "N"]:
    #             tx_intervals.append(
    #                 (-tx_pos, -tx_pos)
    #             )  # negative values to indicate skipped regions
    #             gx_intervals.append((gx_pos, gx_pos + length - 1))
    #             gx_pos += length
    #         elif operation in ["I", "S"]:
    #             tx_intervals.append((tx_pos, tx_pos + length - 1))
    #             gx_intervals.append(
    #                 (-gx_pos, -gx_pos)
    #             )  # negative values to indicate skipped regions
    #             tx_pos += length
    #         else:
    #             pass  # I don't want to deal with H, P

    #     print(f"Genomic intervals: {gx_intervals}")
    #     print(f"Transcript intervals: {tx_intervals}")
    #     # TODO: explore data structures to store intervals
    #     self.tx_intervals = tx_intervals
    #     self.gx_intervals = gx_intervals

    # def get_genomic_coordinates_from_precomputed_intervals(self, tx_start):
    #     """
    #     Given a transcript coordinate (tx_start), return the corresponding genomic coordinate using precomputed intervals.
    #     I don't use this because I don't consider it as effcient but it helps me check my logic.
    #     """
    #     for i, tx_interval in enumerate(self.tx_intervals):
    #         tx_interval_start = tx_interval[0]
    #         tx_interval_end = tx_interval[1]
    #         if tx_start >= tx_interval_start and tx_start <= tx_interval_end:
    #             tx_pos = tx_start - tx_interval_start
    #             genomic_pos = self.gx_intervals[i][0] + tx_pos
    #             return genomic_pos

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
            elif operation in ["D", "N"]:
                gx_pos += length
            elif operation in ["I", "S"]:
                tx_pos += length
            else:
                pass # I don't want to deal with H, P because they don't affect the genomic position
        self.exonIntervals = exonIntervals
        self.tx_len = tx_pos

    def get_genomic_coordinates_from_precomputed_intervals2(self, tx_start: int):
        """
        Given a transcript coordinate (tx_start), return the corresponding genomic coordinate using precomputed intervals.
        """
        if self.exonIntervals is None or self.tx_len == 0:
            raise ValueError("Exon intervals have not been precomputed.")
        if tx_start < 0 or tx_start >= self.tx_len:
            raise ValueError("tx_start {} is out of bounds.".format(tx_start))
        for interval in self.exonIntervals[tx_start]:
            return interval.data + (tx_start - interval.begin)
        return 0  # Return 0 if no interval is found

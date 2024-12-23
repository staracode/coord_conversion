"""
Goal: given a list of transcripts and a list of queries, find the genomic coordinates of the queries in the transcripts.
"""

import sys, os
sys.path.append("./")
from transcript import Transcript
import pandas as pd


def main():
    # Read the data from the file
    transcripts = pd.read_csv(
        "data/transcripts.tsv",
        sep="\t",
        names=["transcript_id", "chromosome", "genomic_start", "cigar"],
    )
    queries = pd.read_csv("data/queries.tsv", sep="\t", names=["query_id", "tx_start"])
    print(transcripts)
    print(queries)


if __name__ == "__main__":
    main()

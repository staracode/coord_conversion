"""
Goal: given a list of transcripts and a list of queries, find the genomic coordinates of the queries in the transcripts.
"""

import sys, os, re

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
    for _, row in transcripts.iterrows():
        transcript = Transcript(
            row["transcript_id"], row["chromosome"], row["genomic_start"], row["cigar"]
        )
        if not transcript.valid_cigar():
            print(
                f"Invalid CIGAR string {transcript.cigar} for transcript {transcript.transcript_id}"
            )
        print(f"Transcript : {transcript.transcript_id}")
        print(f"CIGAR : {transcript.cigar}")
        transcript.visualize_cigar()
        for _, query in queries.iterrows():
            if query["query_id"] != transcript.transcript_id:
                continue
            genomic_start = transcript.get_genomic_coordinate(int(query["tx_start"]))
            print(
                f"{query['query_id']}\t{query['tx_start']}\t{transcript.chromosome}\t{genomic_start}"
            )


if __name__ == "__main__":
    main()

"""
Goal: given a list of transcripts and a list of queries, find the genomic coordinates of the queries in the transcripts.
"""

import sys

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

    # Store each transcript in a dictionary for faster lookup
    all_transcripts = {}
    for _, row in transcripts.iterrows():
        transcript = Transcript(
            row["transcript_id"], row["chromosome"], row["genomic_start"], row["cigar"]
        )
        print(f"Transcript : {transcript.transcript_id}")
        print(f"CIGAR : {transcript.cigar}")
        transcript.visualize_cigar()
        transcript.precompute_tx_intervals()
        transcript.precomputeIntervals2()
        all_transcripts[transcript.transcript_id] = transcript

    # Find the genomic coordinates of the queries
    output = open("output.tsv", "w")
    for _, query in queries.iterrows():
        if query["query_id"] not in all_transcripts:
            print(f"Transcript {query['query_id']} not found")
            continue
        transcript = all_transcripts[query["query_id"]]
        try:
            genomic_start_precomputed = (
                transcript.get_genomic_coordinates_from_precomputed_intervals(
                    int(query["tx_start"])
                )
            )

        except Exception as e:
            print(f"Error: {e}")
            continue
        try:
            genomic_start_precomputed2 = (
                transcript.get_genomic_coordinates_from_precomputed_intervals2(
                    int(query["tx_start"])
                )
            )
        except Exception as e:
            print(f"Error: {e}")
            continue
        print(
            f"A: {query['query_id']}\t{query['tx_start']}\t{transcript.chromosome}\t{genomic_start_precomputed}"
        )
        print(
            f"B: {query['query_id']}\t{query['tx_start']}\t{transcript.chromosome}\t{genomic_start_precomputed2}"
        )
        output.write(
            f"{query['query_id']}\t{query['tx_start']}\t{transcript.chromosome}\t{genomic_start_precomputed2}\n"
        )
    output.close()


if __name__ == "__main__":
    main()

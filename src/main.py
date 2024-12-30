"""
Goal: given a list of transcripts and a list of queries, find the genomic coordinates of the queries in the transcripts.
"""

import sys

sys.path.append("./")
from transcript import Transcript
import pandas as pd
import argparse
import logging


def main():
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename="example.log", encoding="utf-8", level=logging.DEBUG)
    argparser = argparse.ArgumentParser(
        description="Find genomic coordinates of queries in transcripts"
    )
    transcript_file = argparser.add_argument(
        "-t",
        "--transcript_file",
        type=str,
        help="Path to the transcript file",
        default="data/transcripts.tsv",
    )
    query_file = argparser.add_argument(
        "-q",
        "--query_file",
        type=str,
        help="Path to the query file",
        default="data/queries.tsv",
    )
    output_file = argparser.add_argument(
        "-o",
        "--output_file",
        type=str,
        help="Path to the output file",
        default="output/output.tsv",
    )
    debug = argparser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Print debug information",
        default=False,
    )
    args = argparser.parse_args()
    # Read the data from the file
    try:
        transcripts = pd.read_csv(
            args.transcript_file,
            sep="\t",
            names=["transcript_id", "chromosome", "genomic_start", "cigar"],
            dtype={
                "transcript_id": str,
                "chromosome": str,
                "genomic_start": int,
                "cigar": str,
            },
            header=None,
        )
        if transcripts.shape[1] != 4:
            raise ValueError(
                "Transcript file must contain exactly four columns: transcript_id, chromosome, genomic_start, cigar"
            )
        if transcripts.isnull().values.any():
            raise ValueError("Transcript file contains empty values")
        if transcripts["genomic_start"].dtype != int:
            raise ValueError("genomic_start must be an integer")
    except Exception as e:
        logger.error(f"Error reading transcript file {args.transcript_file}: {e}")
        sys.exit(1)

    try:
        queries = pd.read_csv(
            args.query_file,
            sep="\t",
            names=["query_id", "tx_start"],
            header=None,
            dtype={"query_id": str, "tx_start": int},
        )
        if queries.shape[1] != 2:
            raise ValueError(
                "Query file must contain exactly two columns: query_id and tx_start"
            )
        if queries.isnull().values.any():
            raise ValueError("Query file contains empty values")
        if queries["tx_start"].dtype != int:
            raise ValueError("tx_start must be an integer")
    except Exception as e:
        logger.error(f"Error reading query file {args.query_file}: {e}")
        sys.exit(1)

    # Store each transcript in a dictionary for faster lookup
    all_transcripts = {}
    for _, row in transcripts.iterrows():
        try:
            transcript = Transcript(
                row["transcript_id"],
                row["chromosome"],
                row["genomic_start"],
                row["cigar"],
            )
            transcript.precomputeIntervals2()
            all_transcripts[transcript.transcript_id] = transcript
        except Exception as e:
            logger.error(f"Error creating transcript {row['transcript_id']}: {e}")
            sys.exit(1)

    # Process queries
    output = open(args.output_file, "w")
    for _, query in queries.iterrows():
        if query["query_id"] not in all_transcripts:
            logger.warning(f"Transcript {query['query_id']} not found")
            continue
        transcript = all_transcripts[query["query_id"]]

        try:
            genomic_start_precomputed2 = (
                transcript.get_genomic_coordinates_from_precomputed_intervals2(
                    int(query["tx_start"])
                )
            )
            if genomic_start_precomputed2 == 0:
                logger.error(
                    f"Error converting transcript {transcript.transcript_id} to genomic coordinates"
                )
        except Exception as e:
            logger.error(
                f"Error converting transcript {transcript.transcript_id} to genomic coordinates: {e}"
            )
            continue
        output.write(
            f"{query['query_id']}\t{query['tx_start']}\t{transcript.chromosome}\t{genomic_start_precomputed2}\n"
        )
    output.close()


if __name__ == "__main__":
    main()

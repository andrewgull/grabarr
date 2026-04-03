#!/usr/bin/env python3

"""a script for genomes and repeats generation"""

import random
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def get_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generates a random genome with one or more repeat types inserted at\n"
            "specified positions.  Use -r once per repeat type:\n\n"
            "  -r LENGTH:POS1,POS2,...[:DIFF]\n\n"
            "DIFF (optional, 0.0–1.0) is the fraction of positions that differ\n"
            "between copies.  Each copy is mutated independently.\n\n"
            "Examples:\n"
            "  # identical copies\n"
            "  gr_generator.py -gl 100000 -r 50:1000,5000 -o out.fasta\n"
            "  # two types, second has 10% variation between copies\n"
            "  gr_generator.py -gl 100000 -r 50:1000,5000 -r 30:20000,30000:0.1 -o out.fasta"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-gl", "--genome_length", type=int, metavar="<integer>",
                        help="length of the background genome to generate")
    parser.add_argument("-r", "--repeat", type=str, metavar="<length:pos1,pos2[:diff]>",
                        action="append", dest="repeats",
                        help="repeat type: length, positions, optional diff fraction (repeatable)")
    parser.add_argument("-o", "--output", type=str, metavar="<output file>",
                        help="output file, fasta format")
    return parser.parse_args()


def random_seq(length):
    return ''.join(random.choice('ACTG') for _ in range(length))


def mutate_seq(seq, diff_percent):
    """Return a copy of seq with exactly round(len(seq) * diff_percent) positions
    changed to a randomly chosen different nucleotide."""
    n_mut = round(len(seq) * diff_percent)
    if n_mut == 0:
        return seq
    positions = random.sample(range(len(seq)), n_mut)
    bases = list(seq)
    for i in positions:
        bases[i] = random.choice([b for b in "ACTG" if b != bases[i]])
    return "".join(bases)


def insert_repeats(dna, insertions):
    """Insert sequences into dna so each appears at exactly its stated position
    in the final string.

    insertions: list of (position, sequence) tuples, in any order.
    Positions must be non-overlapping in the final string.
    """
    sorted_ins = sorted(insertions, key=lambda x: x[0])
    result = []
    prev_orig = 0
    offset = 0  # cumulative bp inserted before current position
    for pos, seq in sorted_ins:
        orig_pos = pos - offset
        result.append(dna[prev_orig:orig_pos])
        result.append(seq)
        prev_orig = orig_pos
        offset += len(seq)
    result.append(dna[prev_orig:])
    return "".join(result)


def main():
    args = get_args()

    repeat_units = []   # canonical (unmutated) sequence per type
    all_insertions = []

    for spec in (args.repeats or []):
        parts = spec.split(":")
        length = int(parts[0])
        positions = [int(p) for p in parts[1].split(",")]
        diff = float(parts[2]) if len(parts) == 3 else 0.0

        base_seq = random_seq(length)
        repeat_units.append(base_seq)
        for pos in positions:
            all_insertions.append((pos, mutate_seq(base_seq, diff)))

    genome = insert_repeats(random_seq(args.genome_length), all_insertions)

    SeqIO.write([SeqRecord(Seq(genome), id="genomic_string", description="")],
                args.output, "fasta")

    out = Path(args.output)
    repeats_path = out.with_name(out.stem + "_repeats" + out.suffix)
    repeat_records = [
        SeqRecord(Seq(seq), id=f"repeat_unit_{i}",
                  description=f"repeat_unit_{i}_length_{len(seq)}")
        for i, seq in enumerate(repeat_units, start=1)
    ]
    SeqIO.write(repeat_records, repeats_path, "fasta")

    total = len(all_insertions)
    types = len(repeat_units)
    print(f"Generated genome with {total} repeat insertion(s) across {types} type(s)")
    print(f"Repeat units written to {repeats_path}")


if __name__ == "__main__":
    main()

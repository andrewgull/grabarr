import pytest
from Bio import SeqIO
from gr_generator import random_seq, insert_repeats, mutate_seq


# --- random_seq ---

def test_random_seq_returns_correct_length():
    assert len(random_seq(100)) == 100


def test_random_seq_only_valid_nucleotides():
    seq = random_seq(1000)
    assert set(seq) <= set("ACTG")


def test_random_seq_zero_length():
    assert random_seq(0) == ""


# --- insert_repeats ---
# insert_repeats(dna, insertions) where insertions = [(position, sequence), ...]

def test_insert_repeats_single_insertion_at_exact_position():
    dna = "A" * 100
    repeat = "TTTT"
    result = insert_repeats(dna, [(10, repeat)])
    assert result[10:14] == repeat


def test_insert_repeats_same_type_multiple_positions_are_exact():
    dna = "A" * 200
    repeat = "GCGC"
    positions = [10, 50, 100]
    result = insert_repeats(dna, [(p, repeat) for p in positions])
    for pos in positions:
        assert result[pos:pos + len(repeat)] == repeat, \
            f"Repeat not found at position {pos}"


def test_insert_repeats_genome_grows_by_sum_of_all_repeat_lengths():
    original_len = 200
    insertions = [(10, "GGGG"), (50, "GGGG"), (100, "GGGG")]
    result = insert_repeats("A" * original_len, insertions)
    assert len(result) == original_len + sum(len(seq) for _, seq in insertions)


def test_insert_repeats_flanks_unchanged():
    dna = "A" * 100
    result = insert_repeats(dna, [(5, "TT")])
    assert result[:5] == "AAAAA"
    assert result[7:] == "A" * 95


def test_insert_repeats_empty_insertions_returns_unchanged():
    dna = "ACGT" * 10
    assert insert_repeats(dna, []) == dna


def test_insert_repeats_two_different_types_at_exact_positions():
    """Two repeat types of different length both land at their stated positions."""
    dna = "A" * 300
    insertions = [(50, "TTTT"), (150, "GCGCGCGC")]  # lengths 4 and 8
    result = insert_repeats(dna, insertions)
    assert result[50:54] == "TTTT"
    assert result[150:158] == "GCGCGCGC"


def test_insert_repeats_different_type_positions_exact_unsorted():
    """Insertions given in reverse order should still land at correct positions."""
    dna = "A" * 300
    insertions = [(150, "GCGC"), (50, "TTTT")]  # deliberately reversed
    result = insert_repeats(dna, insertions)
    assert result[50:54] == "TTTT"
    assert result[150:154] == "GCGC"


def test_insert_repeats_genome_grows_by_sum_of_different_length_repeats():
    dna = "A" * 500
    insertions = [(100, "TT"), (200, "GCGCGC"), (300, "AAAA")]
    result = insert_repeats(dna, insertions)
    assert len(result) == 500 + 2 + 6 + 4


# --- mutate_seq ---

def test_mutate_seq_zero_diff_returns_identical_sequence():
    seq = "ACGTACGTACGT"
    assert mutate_seq(seq, 0.0) == seq


def test_mutate_seq_changes_exact_number_of_positions():
    seq = "A" * 100
    result = mutate_seq(seq, 0.2)
    n_changed = sum(a != b for a, b in zip(seq, result))
    assert n_changed == round(100 * 0.2)


def test_mutate_seq_mutated_positions_have_different_nucleotide():
    seq = "AAAAAAAAAA"
    result = mutate_seq(seq, 0.3)
    for orig, mut in zip(seq, result):
        if orig != mut:
            assert mut in "ACTG"
            assert mut != orig


def test_mutate_seq_unmutated_positions_match_original():
    seq = "ACGTACGTACGTACGT"
    result = mutate_seq(seq, 0.25)
    for orig, mut in zip(seq, result):
        if orig == mut:
            # position unchanged — fine
            pass
    # count unchanged positions
    n_unchanged = sum(a == b for a, b in zip(seq, result))
    assert n_unchanged == len(seq) - round(len(seq) * 0.25)


def test_mutate_seq_output_length_unchanged():
    seq = random_seq(200)
    assert len(mutate_seq(seq, 0.15)) == 200


def test_mutate_seq_only_valid_nucleotides():
    seq = random_seq(100)
    result = mutate_seq(seq, 0.1)
    assert set(result) <= set("ACTG")


# --- integration: main writes a valid FASTA ---

def test_main_single_type_writes_genome_and_repeat_records(tmp_path):
    import sys
    outfile = tmp_path / "genome.fasta"
    sys.argv = ["gr_generator.py", "-gl", "500",
                "-r", "10:50,100,200",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    assert records[0].id == "genomic_string"
    assert records[1].id == "repeat_unit_1"
    assert len(records[1].seq) == 10


def test_main_single_type_repeat_appears_at_input_positions(tmp_path):
    import sys
    outfile = tmp_path / "genome.fasta"
    positions = [50, 150, 300]
    sys.argv = ["gr_generator.py", "-gl", "1000",
                "-r", f"20:{','.join(str(p) for p in positions)}",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    genome = str(records[0].seq)
    repeat = str(records[1].seq)
    for pos in positions:
        assert genome[pos:pos + len(repeat)] == repeat, \
            f"Repeat not found at input position {pos}"


def test_main_multi_type_writes_genome_plus_one_record_per_type(tmp_path):
    import sys
    outfile = tmp_path / "genome.fasta"
    sys.argv = ["gr_generator.py", "-gl", "1000",
                "-r", "20:100,200",
                "-r", "40:400,600",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    # 1 genomic_string + 2 repeat_unit records
    assert len(records) == 3
    assert records[0].id == "genomic_string"
    assert records[1].id == "repeat_unit_1"
    assert records[2].id == "repeat_unit_2"


def test_main_multi_type_repeat_lengths_are_correct(tmp_path):
    import sys
    outfile = tmp_path / "genome.fasta"
    sys.argv = ["gr_generator.py", "-gl", "1000",
                "-r", "20:100,200",
                "-r", "40:400,600",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    assert len(records[1].seq) == 20
    assert len(records[2].seq) == 40


def test_main_multi_type_each_repeat_at_exact_position(tmp_path):
    import sys
    outfile = tmp_path / "genome.fasta"
    sys.argv = ["gr_generator.py", "-gl", "2000",
                "-r", "20:100,300",
                "-r", "50:600,900",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    genome = str(records[0].seq)
    repeat1 = str(records[1].seq)
    repeat2 = str(records[2].seq)

    for pos in [100, 300]:
        assert genome[pos:pos + 20] == repeat1, f"repeat1 not at position {pos}"
    for pos in [600, 900]:
        assert genome[pos:pos + 50] == repeat2, f"repeat2 not at position {pos}"


# --- diff (mutation) feature ---

def test_main_zero_diff_all_copies_identical_in_genome(tmp_path):
    import sys
    outfile = tmp_path / "genome.fasta"
    sys.argv = ["gr_generator.py", "-gl", "2000",
                "-r", "40:100,500,1000:0.0",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    genome = str(records[0].seq)
    base = str(records[1].seq)

    for pos in [100, 500, 1000]:
        assert genome[pos:pos + 40] == base


def test_main_nonzero_diff_repeat_unit_record_is_canonical_base(tmp_path):
    """The repeat_unit record always holds the unmutated base sequence."""
    import sys
    outfile = tmp_path / "genome.fasta"
    sys.argv = ["gr_generator.py", "-gl", "2000",
                "-r", "50:200,700:0.2",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    assert len(records[1].seq) == 50


def test_main_nonzero_diff_copies_differ_from_base_sequence(tmp_path):
    """With diff>0, inserted copies should differ from the canonical base."""
    import sys
    outfile = tmp_path / "genome.fasta"
    sys.argv = ["gr_generator.py", "-gl", "2000",
                "-r", "100:200,700:0.2",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    genome = str(records[0].seq)
    base = str(records[1].seq)

    copies = [genome[200:300], genome[700:800]]
    n_diffs = [sum(a != b for a, b in zip(base, copy)) for copy in copies]
    # each copy should have exactly round(100 * 0.2) = 20 positions changed
    assert all(d == round(100 * 0.2) for d in n_diffs)


def test_main_nonzero_diff_copies_at_exact_positions(tmp_path):
    """Mutated copies still start at the stated positions."""
    import sys
    outfile = tmp_path / "genome.fasta"
    sys.argv = ["gr_generator.py", "-gl", "2000",
                "-r", "30:100,400:0.1",
                "-o", str(outfile)]

    from gr_generator import main
    main()

    records = list(SeqIO.parse(str(outfile), "fasta"))
    genome = str(records[0].seq)
    base = str(records[1].seq)

    for pos in [100, 400]:
        copy = genome[pos:pos + 30]
        n_diff = sum(a != b for a, b in zip(base, copy))
        # exactly round(30 * 0.1) = 3 positions should differ
        assert n_diff == round(30 * 0.1), \
            f"copy at {pos} has {n_diff} diffs, expected {round(30 * 0.1)}"

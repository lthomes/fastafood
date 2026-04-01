from typing import Iterable
from .sequence import DNASequence


def to_fasta(sequences: Iterable[DNASequence], file_path: str, line_width: int = 60):
    with open(file_path, "w") as f:
        for i, seq in enumerate(sequences):
            header = seq.name if seq.name else f"seq_{i}"
            f.write(f">{header}\n")

            s = seq.sequence
            for j in range(0, len(s), line_width):
                f.write(s[j:j+line_width] + "\n")

def read_fasta(file_path: str):
    sequences = []

    with open(file_path, "r") as f:
        name = None
        seq_chunks = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if name is not None:
                    sequences.append(
                        DNASequence("".join(seq_chunks), name=name)
                    )
                name = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        # dernière séquence
        if name is not None:
            sequences.append(
                DNASequence("".join(seq_chunks), name=name)
            )

    return sequences

def count_fasta(file_path: str) -> int:
    count = 0
    with open(file_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

def list_fasta_headers(file_path: str):
    headers = []
    with open(file_path) as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line.strip()[1:])
    return headers

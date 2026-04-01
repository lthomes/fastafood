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
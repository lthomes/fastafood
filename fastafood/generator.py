from typing import Iterable, Set
from .sequence import DNASequence

class VariantGenerator:

    def __init__(self, dna: DNASequence):
        self.dna = dna

    def generate_snps(self, protected_positions: Set[int] = None) -> Iterable[DNASequence]:
        if protected_positions is None:
            protected_positions = set()

        for i, base in enumerate(self.dna.sequence):
            if i in protected_positions:
                continue

            for alt in DNASequence.VALID_BASES:
                if alt != base:
                    yield self.dna.mutate(i, alt)

    def generate_deletions(self, sizes=(1,), protected_positions: Set[int] = None, step: int = None):
        if protected_positions is None:
            protected_positions = set()

        seq = self.dna.sequence
        n = len(seq)

        for size in sizes:
            # Si step est None, on glisse de 1. Sinon, on utilise la valeur fournie.
            actual_step = step if step is not None else 1
            
            for i in range(0, n - size + 1, actual_step):
                if any(pos in protected_positions for pos in range(i, i + size)):
                    continue

                new_seq = seq[:i] + seq[i + size:]
                yield DNASequence(new_seq)

    def generate_duplications(self, sizes=(1,), times=2, protected_positions: Set[int] = None):
        if protected_positions is None:
            protected_positions = set()

        seq = self.dna.sequence
        n = len(seq)

        for size in sizes:
            for i in range(n - size + 1):
                if any(pos in protected_positions for pos in range(i, i + size)):
                    continue

                fragment = seq[i:i+size]
                duplicated = fragment * times
                new_seq = seq[:i] + duplicated + seq[i+size:]
                yield DNASequence(new_seq)
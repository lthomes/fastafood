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

    def generate_insertion(self, pos: int, insert_seq: str) -> DNASequence:
        # On vérifie si la position est dans la séquence
        if pos < 0 or pos > len(self.dna.sequence):
            raise IndexError("Position d'insertion hors limites")
            
        return self.dna.insert(pos, insert_seq)

    def generate_multiple_insertions(self, insertions: list) -> DNASequence:
        """
        insertions: liste de tuples (pos, seq_to_ins)
        """
        # On trie par position décroissante pour ne pas décaler les index
        # pendant qu'on reconstruit la chaîne.
        sorted_insertions = sorted(insertions, key=lambda x: x[0], reverse=True)
        
        current_seq_str = self.dna.sequence
        for pos, seq_to_ins in sorted_insertions:
            # Validation de la séquence à insérer
            DNASequence(seq_to_ins) 
            # Insertion directe par slicing
            current_seq_str = current_seq_str[:pos] + seq_to_ins.upper() + current_seq_str[pos:]
            
        return DNASequence(current_seq_str, name=f"{self.dna.name}_multi_ins")

    def generate_inversions(self, sizes=(2,), step: int = None, protected_positions: Set[int] = None):
        if protected_positions is None:
            protected_positions = set()

        seq = self.dna.sequence
        n = len(seq)

        for size in sizes:
            actual_step = step if step is not None else 1
            
            for i in range(0, n - size + 1, actual_step):
                # On vérifie si une partie de la région à inverser est protégée
                if any(pos in protected_positions for pos in range(i, i + size)):
                    continue

                # On extrait, on inverse, et on réinsère
                fragment = seq[i:i+size]
                inverted_fragment = fragment[::-1]
                new_seq = seq[:i] + inverted_fragment + seq[i+size:] 
                
                yield DNASequence(new_seq, name=f"{self.dna.name}_inv_{i}_{size}")
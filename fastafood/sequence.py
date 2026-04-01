
from typing import List

class DNASequence:
    VALID_BASES = {"A", "T", "C", "G"}

    COMPLEMENT = str.maketrans("ATCG", "TAGC")

    CODON_TABLE = {
        "ATA":"I","ATC":"I","ATT":"I","ATG":"M",
        "ACA":"T","ACC":"T","ACG":"T","ACT":"T",
        "AAC":"N","AAT":"N","AAA":"K","AAG":"K",
        "AGC":"S","AGT":"S","AGA":"R","AGG":"R",

        "CTA":"L","CTC":"L","CTG":"L","CTT":"L",
        "CCA":"P","CCC":"P","CCG":"P","CCT":"P",
        "CAC":"H","CAT":"H","CAA":"Q","CAG":"Q",
        "CGA":"R","CGC":"R","CGG":"R","CGT":"R",

        "GTA":"V","GTC":"V","GTG":"V","GTT":"V",
        "GCA":"A","GCC":"A","GCG":"A","GCT":"A",
        "GAC":"D","GAT":"D","GAA":"E","GAG":"E",
        "GGA":"G","GGC":"G","GGG":"G","GGT":"G",

        "TCA":"S","TCC":"S","TCG":"S","TCT":"S",
        "TTC":"F","TTT":"F","TTA":"L","TTG":"L",
        "TAC":"Y","TAT":"Y","TAA":"*","TAG":"*",
        "TGC":"C","TGT":"C","TGA":"*","TGG":"W",
    }

    def __init__(self, sequence: str, name: str = None):
        self.sequence = sequence.upper()
        self.name = name
        self._validate()

    def _validate(self):
        if not set(self.sequence).issubset(self.VALID_BASES):
            raise ValueError("Invalid DNA sequence")

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return f"DNASequence('{self.sequence}')"

    def __eq__(self, other):
        return isinstance(other, DNASequence) and self.sequence == other.sequence

    # --- Mutations ---
    def mutate(self, pos: int, new_base: str):
        if new_base not in self.VALID_BASES:
            raise ValueError("Invalid base")
        seq = list(self.sequence)
        seq[pos] = new_base
        return DNASequence("".join(seq), self.name)

    # --- Trimming ---
    def trim_5prime(self, n: int):
        return DNASequence(self.sequence[n:], self.name)

    def trim_3prime(self, n: int):
        return DNASequence(self.sequence[:-n] if n > 0 else self.sequence, self.name)

    def trim(self, n5: int = 0, n3: int = 0):
        return DNASequence(self.sequence[n5: len(self.sequence) - n3], self.name)

    # --- Reverse complement ---
    def complement(self):
        return DNASequence(self.sequence.translate(self.COMPLEMENT), self.name)

    def reverse(self):
        return DNASequence(self.sequence[::-1], self.name)

    def reverse_complement(self):
        return DNASequence(self.sequence.translate(self.COMPLEMENT)[::-1], self.name)

    # --- Translation ---
    def translate(self, frame: int = 0, stop_at_stop: bool = False) -> str:
        seq = self.sequence[frame:]
        protein: List[str] = []

        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa = self.CODON_TABLE.get(codon, "X")

            if stop_at_stop and aa == "*":
                break

            protein.append(aa)

        return "".join(protein)
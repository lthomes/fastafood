          .----------------.
      _.-'    FASTA FOOD    '-._
    .'__________________________'.
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      |   _  ..  _  ..  _  ..  |
      |  / \/  \/ \/  \/ \/  \ |
      |  \ /\  /\ /\  /\ /\  / |
      |   '  ''  '  ''  '  ''  |
    !____________________________!
    '----------------------------'

        --- fastafood v0.1.0 ---
         "From PASTA to FASTA."

## Features

- SNP generation (with protected positions)
- Deletions (base, codon, etc.)
- Duplications
- Reverse complement
- Translation
- FASTA export

## Example

```python
from fastafood import DNASequence, VariantGenerator, to_fasta

dna = DNASequence("ATGC")

gen = VariantGenerator(dna)

variants = list(gen.generate_snps(protected_positions={0}))

to_fasta(variants, "variants.fasta")


---

## 🧪 4. Installer en local (mode dev)

Dans le dossier racine :

```bash
pip install -e .
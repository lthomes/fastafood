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
      --- fastafood v0.2.0 ---
       "From PASTA to FASTA."

**fastafood** est un couteau suisse bio-informatique en Python. Il permet de manipuler des séquences ADN, de simuler des mutations complexes et de générer des bibliothèques de variants avec une précision chirurgicale.

## Installation

```shell
git clone [https://github.com/lthomes/fastafood.git](https://github.com/lthomes/fastafood.git)
cd fastafood
pip install .
```

## Utilisation (CLI)

### 1. Entrée des données

L'outil accepte des fichiers FASTA ou des chaînes de caractères brutes :

Fichier : `fastafood --file genome.fasta`

Séquence brute : `fastafood --seq ATGCATGC`

### 2. Transformations de base

Reverse Complement : `--revcomp`

Trimming : `--trim 5 10` (supprime 5 bases en 5' et 10 bases en 3')

Traduction Protéique : `--translate` ou `--scan-translate` (3 cadres).

Gestion du Stop : `--stop yes` (tronque au premier \*).

## Génération de Variants

### A. Mutations Automatisées (Massives)

Générez tous les variants possibles pour une taille donnée :

SNPs : `--snps`

Délétions : `--deletions 2 3` (génère toutes les délétions de taille 2 et 3)

Duplications : `--duplications 5`

Inversions : `--inversions 10`

Options : `--step 3` (définit le saut de la fenêtre glissante).

### B. Mutations Ciblées (Précision)

Effectuez des modifications spécifiques à des positions précises :

SNP ciblé : `--snps-target 6 T` (remplace la base 6 par T)

Remplacement de zone : `--snps-target 6 10 TTATT` (remplace du nucléotide 6 à 10 par TTATT)

Délétion ciblée : `--deletions-target 2 5` (supprime les bases 2 à 5)

Insertion : `--insertion 4 --ins-seq ATGC` (insère ATGC à la position 4)

### C. Duplications et Réarrangements

Duplication simple : `--duplications-target 2 5` (duplique la zone 2-5 juste après elle-même)

Duplication distante : `--duplications-target 2 5 12` (duplique 2-5 et l'insère à la position 12)

Options avancées :

- `--reversed yes` : inverse le fragment avant insertion.

- `--times 3` : insère le fragment en 3 exemplaires.

## Système de Protection

Rendez certaines zones "intouchables" pour les générateurs de mutations :

Ponctuel : `--protect 4 8 15` (protège les bases spécifiées)

Plages : `--protect-range 10 20 50 60` (protège de 10 à 20 et de 50 à 60)

Extrémités

- `--protect-extremity 15` (protège les 15 premières bases)

- `--protect-extremity 10 80` (protège les 10 premières bases ET de la base 80 jusqu'à la fin)

## Export & Formatage

Sortie fichier : `--out results.fasta`

Formatage : `--format split` (découpage à 60 car.) ou `one-line` (par défaut).

## Utilisation comme Bibliothèque

Intégrez la logique de fastafood dans vos propres scripts :

```Python
from fastafood import DNASequence, VariantGenerator

# Initialisation
dna = DNASequence("ATGCGTACGTAG", name="wild_type")

# Protection des 3 premiers nucléotides
protected = {0, 1, 2}

# Génération ciblée
gen = VariantGenerator(dna)
variant = gen.generate_targeted_duplication(start=0, end=4, target=10, times=2)

print(f"Original : {dna.sequence}")
print(f"Variant  : {variant.sequence}")
```

## Points clés

**1-based Indexing** : Les positions en CLI correspondent aux standards bio-informatiques (la première base est la n°1).

**Validation stricte** : Refuse les caractères non-ADN lors des insertions/remplacements.

**Auto-naming** : Génère des noms explicites pour les variants (seq_dup_x3, seq_rev_comp, etc.).

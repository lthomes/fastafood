import argparse
import sys
import os
from types import SimpleNamespace
from fastafood import read_fasta, to_fasta, VariantGenerator, DNASequence

# "Magic" line pour ANSI sur Windows
os.system('') 

def main():
    parser = argparse.ArgumentParser(description="Fastafood CLI - Mode Fichier ou Séquence brute")
    if len(sys.argv) == 1 or "--help" in sys.argv or "-h" in sys.argv:
        fast_a_food_splash()
        if len(sys.argv) == 1: return
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--file", help="Fichier FASTA d'entrée")
    input_group.add_argument("--seq", help="Séquence ADN brute (ex: ATGC...)")

    parser.add_argument("--out", help="Fichier de sortie (optionnel)")
    
    # --- Options de Formatage ---
    fmt_group = parser.add_argument_group("Formatage")
    fmt_group.add_argument("--format", choices=["one-line", "split"], default="one-line", 
                           help="Format de la séquence (défaut: one-line)")

    # --- Transformations & Mutations ---
    transform = parser.add_argument_group("Transformations")
    transform.add_argument("--revcomp", action="store_true", help="Reverse complement")
    transform.add_argument("--trim", nargs=2, type=int, metavar=('5p', '3p'))
    transform.add_argument("--translate", action="store_true", help="Traduit en protéine")
    transform.add_argument("--scan-translate", action="store_true", help="Scanning translation (6 cadres)")
    # Nouvel argument pour la gestion du STOP
    transform.add_argument("--stop", choices=["yes", "no"], default="no",
                           help="Arrêter la traduction au premier STOP rencontré (défaut: no)")

    mut = parser.add_argument_group("Mutations")
    mut.add_argument("--snps", action="store_true")
    mut.add_argument("--deletions", nargs="+", type=int, metavar="SIZE")
    mut.add_argument("--step", type=int)
    mut.add_argument("--duplications", nargs="+", type=int, metavar="SIZE")
    mut.add_argument("--protect", nargs="+", type=int)

    args = parser.parse_args()

    # Conversion de l'argument --stop en booléen pour les méthodes
    stop_flag = True if args.stop == "yes" else False

    # 1. Acquisition
    input_sequences = []
    if args.file:
        try:
            input_sequences = read_fasta(args.file)
        except Exception as e:
            print(f"Erreur fichier : {e}")
            sys.exit(1)
    elif args.seq:
        input_sequences = [DNASequence(args.seq, name="manual_input")]

    final_results = []
    one_based = lambda pos: pos - 1
    protected = set(one_based(pos) for pos in args.protect) if args.protect else set() 

    for seq in input_sequences:
        # 2. Transformations ADN
        if args.revcomp:
            seq = seq.reverse_complement()
        if args.trim:
            seq = seq.trim(n5=args.trim[0], n3=args.trim[1]) 

        # 3. Variants
        gen = VariantGenerator(seq)
        variants = [seq]

        if args.snps:
            variants.extend(list(gen.generate_snps(protected_positions=protected))) 
        if args.deletions:
            variants.extend(list(gen.generate_deletions(sizes=args.deletions, protected_positions=protected, step=args.step)))
        if args.duplications:
            variants.extend(list(gen.generate_duplications(sizes=args.duplications, protected_positions=protected)))

        # 4. Traduction Exclusive
        def get_safe_name(v):
            return v.name if v.name else "sequence"

        if args.translate:
            for v in variants:
                prot_seq = v.translate(stop_at_stop=stop_flag)
                final_results.append({"name": f"{get_safe_name(v)}_prot", "sequence": prot_seq})
        elif args.scan_translate:
            for v in variants:
                prot_seqs = v.scanning_translation(stop_at_stop=stop_flag)
                for i, prot in enumerate(prot_seqs, 1):
                    final_results.append({"name": f"{get_safe_name(v)}_frame{i}_scan", "sequence": prot})
        else:
            final_results.extend(variants)

    # 5. Export / Affichage
    width = 60 if args.format == "split" else 10**9 
    
    export_list = []
    for index, item in enumerate(final_results):
        if isinstance(item, DNASequence):
            raw_name = item.name
            sequence = item.sequence
        else:
            raw_name = item.get("name")
            sequence = item.get("sequence")
        
        if not raw_name:
            raw_name = "sequence"

        indexed_name = f"{raw_name.lstrip('>')}_{index}"
        obj = SimpleNamespace(name=indexed_name, sequence=sequence)
        export_list.append(obj)

    if args.out:
        to_fasta(export_list, args.out, line_width=width) 
        print(f"Succès : {len(export_list)} séquences ({args.format}) écrites dans {args.out}")
    else:
        for res in export_list:
            print(f">{res.name}")
            s = res.sequence
            if args.format == "split":
                for i in range(0, len(s), 60):
                    print(s[i:i+60])
            else:
                print(s)

def fast_a_food_splash():
    BUN, LETTUCE, DNA = "\033[38;5;214m", "\033[32m", "\033[36m"
    RESET, BOLD, KETCHUP = "\033[0m", "\033[1m", "\033[31m"
    print(f"""
          {BUN}.----------------.{RESET}
      {BUN}_.-'    {BOLD}{KETCHUP}FASTA FOOD{RESET}{BUN}    '-._{RESET}
    {BUN}.'__________________________'.{RESET}
    {LETTUCE}!~~~~~~~~~~~~~~~~~~~~~~~~~~~~!{RESET}
      {DNA}|   _  ..  _  ..  _  ..  |{RESET}
      {DNA}|  / \/  \/ \/  \/ \/  \ |{RESET}
      {DNA}|  \ /\  /\ /\  /\ /\  / |{RESET}
      {DNA}|   '  ''  '  ''  '  ''  |{RESET}
    {BUN}!____________________________!{RESET}
    {BUN}'----------------------------'{RESET}
        {BOLD}--- fastafood v0.1.0 ---{RESET}
    """)

if __name__ == "__main__":
    main()
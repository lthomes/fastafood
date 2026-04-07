import argparse
import sys
import os
from types import SimpleNamespace
from fastafood import read_fasta, to_fasta, VariantGenerator, DNASequence

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
    transform.add_argument("--stop", choices=["yes", "no"], default="no",
                           help="Arrêter la traduction au premier STOP rencontré (défaut: no)")

    mut = parser.add_argument_group("Mutations")
    mut.add_argument("--snps", action="store_true")
    mut.add_argument("--deletions", nargs="+", type=int, metavar="SIZE")
    mut.add_argument("--step", type=int)
    mut.add_argument("--duplications", nargs="+", type=int, metavar="SIZE")
    mut.add_argument("--protect", nargs="+", type=int)
    mut.add_argument("--inversions", nargs="*", type=int, metavar="SIZE", 
                     help="Génère des inversions. Si aucune taille n'est fournie, utilise 2.")
    mut.add_argument("--insertion", nargs="+", type=int, metavar="POS", 
                     help="Position(s) d'insertion (1-based)")
    mut.add_argument("--ins-seq", nargs="+", type=str, metavar="SEQ", 
                     help="Séquence(s) ADN à insérer")
    mut.add_argument("--snps-target", nargs="+", help="Cible : 'POS BASE' ou 'START END SEQ' (1-based)")
    mut.add_argument("--deletions-target", nargs=2, type=int, metavar=('START', 'END'), help="Supprime de START à END (1-based)")
    mut.add_argument("--duplications-target", nargs="+", type=int, 
                 help="Cible : 'START END' ou 'START END TARGET' (1-based)")
    mut.add_argument("--reversed", choices=["yes", "no"], default="no",
                    help="Inverse la séquence dupliquée avant insertion (défaut: no)")
    mut.add_argument("--times", type=int, default=1, help="Nombre de répétitions pour la duplication (défaut: 1)")
    mut.add_argument("--protect-range", nargs="+", type=int, 
                 help="Liste de START STOP à protéger (ex: 3 5 12 13)")
    mut.add_argument("--protect-extremity", nargs="+", type=int, 
                    help="Protège de 1 à POS, ou 1 à POS1 ET POS2 à FIN")

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

        # On initialise le set de protection pour CETTE séquence
        protected = set()

        # 1. Protection ponctuelle existante (--protect)
        if args.protect:
            protected.update(p - 1 for p in args.protect)

        # 2. Protection par plages (--protect-range)
        if args.protect_range:
            # On itère par bonds de 2 pour attraper les couples (start, stop)
            for i in range(0, len(args.protect_range), 2):
                try:
                    start = args.protect_range[i]
                    stop = args.protect_range[i+1]
                    # On ajoute toutes les positions entre start et stop inclus (1-based)
                    protected.update(range(start - 1, stop))
                except IndexError:
                    print(f"Attention : argument impair pour --protect-range, la dernière valeur a été ignorée.")

        # 3. Protection des extrémités (--protect-extremity)
        if args.protect_extremity:
            if len(args.protect_extremity) == 1:
                # Protège du début jusqu'à POS
                pos = args.protect_extremity[0]
                protected.update(range(0, pos))
            elif len(args.protect_extremity) >= 2:
                # Protège du début à POS1 ET de POS2 à la fin
                pos1 = args.protect_extremity[0]
                pos2 = args.protect_extremity[1]
                protected.update(range(0, pos1))
                protected.update(range(pos2 - 1, len(seq)))

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
        if args.inversions is not None:
            # Si l'utilisateur a mis juste --inversions, on prend [2] par défaut
            inv_sizes = args.inversions if len(args.inversions) > 0 else [2]
            variants.extend(list(gen.generate_inversions(sizes=inv_sizes, protected_positions=protected, step=args.step)))

        # 3. Variants - Gestion des insertions simultanées
        if args.insertion is not None:
            if not args.ins_seq:
                print("Erreur : --insertion nécessite --ins-seq")
                sys.exit(1)

            positions = args.insertion
            sequences = args.ins_seq
            insert_map = []

            # Cas 1 : Une seule séquence pour plusieurs positions
            if len(sequences) == 1:
                insert_map = [(pos - 1, sequences[0]) for pos in positions]
            
            # Cas 2 : Autant de séquences que de positions
            elif len(sequences) == len(positions):
                insert_map = [(p - 1, s) for p, s in zip(positions, sequences)]
            
            # Cas 3 : Déséquilibre
            else:
                print(f"Erreur : Déséquilibre ({len(positions)} pos vs {len(sequences)} seq)")
                sys.exit(1)

            try:
                # On génère UN SEUL variant avec toutes les modifs
                new_variant = gen.generate_multiple_insertions(insert_map)
                variants.append(new_variant)
            except Exception as e:
                print(f"Erreur lors des insertions multiples : {e}")

        # --- Gestion des cibles spécifiques ---

        # SNP ou Remplacement de zone
        if args.snps_target:
            try:
                if len(args.snps_target) == 2:
                    # Format: POS BASE (ex: 6 T)
                    pos = int(args.snps_target[0]) - 1
                    base = args.snps_target[1]

                    if pos in protected:
                        print(f"Erreur : la position {pos + 1} est protégée et ne peut pas être modifiée.")
                    else:
                        variants.append(gen.generate_targeted_snp(pos, pos + 1, base))
                elif len(args.snps_target) == 3:
                    # Format: START END SEQ (ex: 6 10 TTATT)
                    start = int(args.snps_target[0]) - 1
                    end = int(args.snps_target[1])
                    seq_replacement = args.snps_target[2]

                    if pos in protected:
                        print(f"Erreur : la position {pos + 1} est protégée et ne peut pas être modifiée.")
                    else:
                        variants.append(gen.generate_targeted_snp(start, end, seq_replacement))
                    
                else:
                    print("Erreur : --snps-target attend 2 ou 3 arguments.")
            except Exception as e:
                print(f"Erreur --snps-target : {e}")

        # Délétion ciblée
        if args.deletions_target:
            try:
                start = args.deletions_target[0] - 1
                end = args.deletions_target[1]
                if pos in protected:
                    print(f"Erreur : la position {pos + 1} est protégée et ne peut pas être modifiée.")
                else:
                    variants.append(gen.generate_targeted_deletion(start, end))
            except Exception as e:
                print(f"Erreur --deletions-target : {e}")

        # --- Duplication ciblée avec option Reverse ---
        if args.duplications_target:
            try:
                is_reversed = (args.reversed == "yes")
                repeat_count = args.times 
                
                start = args.duplications_target[0] - 1
                end = args.duplications_target[1] 
                
                target_pos = None
                if len(args.duplications_target) == 3:
                    target_pos = args.duplications_target[2] - 1
                    
                new_variant = gen.generate_targeted_duplication(
                    start=start, 
                    end=end, 
                    target=target_pos, 
                    reversed_frag=is_reversed,
                    times=repeat_count
                )
                
                # Nommage dynamique
                tag = "rev_dup" if is_reversed else "dup"
                new_variant.name = f"{seq.name}_{tag}_x{repeat_count}"
                if pos in protected:
                    print(f"Erreur : la position {pos + 1} est protégée et ne peut pas être modifiée.")
                else:
                    variants.append(new_variant)
                
            except Exception as e:
                print(f"Erreur --duplications-target : {e}")

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
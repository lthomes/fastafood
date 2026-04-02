import argparse
import sys
from fastafood import read_fasta, to_fasta, VariantGenerator, DNASequence
import os
os.system('') # This "magic" line enables ANSI escape characters in many Windows shells

def main():
    parser = argparse.ArgumentParser(description="Fastafood CLI - Mode Fichier ou Séquence brute")
    if len(sys.argv) == 1 or "--help" in sys.argv or "-h" in sys.argv:
        fast_a_food_splash()
    
    # Groupe exclusif : Fichier OU Séquence [cite: 13, 14]
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--file", help="Fichier FASTA d'entrée")
    input_group.add_argument("--seq", help="Séquence ADN brute (ex: ATGC...)")

    parser.add_argument("--out", help="Fichier de sortie (optionnel)")
    
    # --- Transformations & Mutations (identique à précédemment) ---
    transform = parser.add_argument_group("Transformations")
    transform.add_argument("--revcomp", action="store_true", help="Reverse complement")
    transform.add_argument("--trim", nargs=2, type=int, metavar=('5p', '3p'))
    transform.add_argument("--translate", action="store_true")

    mut = parser.add_argument_group("Mutations")
    mut.add_argument("--snps", action="store_true", help="Génère tous les SNPs possibles")
    mut.add_argument("--deletions", nargs="+", type=int, metavar="SIZE", help="Génère des délétions")
    mut.add_argument("--step", type=int, help="Pas de glissement pour les délétions/duplications")
    mut.add_argument("--duplications", nargs="+", type=int, metavar="SIZE", help="Génère des duplications")
    mut.add_argument("--protect", nargs="+", type=int, help="Positions à protéger (1-indexed)")

    args = parser.parse_args()

    # 1. Acquisition des séquences
    input_sequences = []
    
    if args.file:
        try:
            input_sequences = read_fasta(args.file) # [cite: 2, 5]
        except Exception as e:
            print(f"Erreur fichier : {e}")
            sys.exit(1)
    elif args.seq:
        # On crée un objet DNASequence à la volée 
        try:
            input_sequences = [DNASequence(args.seq, name="manual_input")]

        except ValueError as e:
            print(f"Erreur séquence : {e}") # Déclenché par _validate 
            sys.exit(1)

    final_results = []
    one_based = lambda pos: pos - 1
    protected = set(one_based(pos) for pos in args.protect) if args.protect else set() 

    for seq in input_sequences:
        # 2. Application des filtres de base
        
        
        if args.revcomp:
            seq = seq.reverse_complement()

        if args.trim:
            seq = seq.trim(n5=args.trim[0], n3=args.trim[1]) 

        # 3. Génération des variants
        gen = VariantGenerator(seq)
        variants = [seq] # On garde l'originale transformée

        if args.snps:
            variants.extend(list(gen.generate_snps(protected_positions=protected))) 
        
        if args.deletions:
        # On passe le paramètre step récupéré de la ligne de commande
            variants.extend(list(gen.generate_deletions(
                sizes=args.deletions, 
                protected_positions=protected, 
                step=args.step
            )))
            
        if args.duplications:
            variants.extend(list(gen.generate_duplications(sizes=args.duplications, protected_positions=protected)))

        # 4. Traduction finale si demandée
        if args.translate:
            for v in variants:
                prot_seq = v.translate(stop_at_stop=True)
                # On crée une pseudo-séquence pour l'export (attention: VALID_BASES bloquerait une DNASequence de AA)
                final_results.append({"name": f"{v.name}_prot", "sequence": prot_seq})
        else:
            final_results.extend(variants)

    # 5. Export ou Affichage
    if args.out:
        # On adapte les objets pour to_fasta
        export_objs = []
        for item in final_results:
            if isinstance(item, DNASequence):
                export_objs.append(item)
            else:
                # Pour les protéines traduites qui sont des dicts
                export_objs.append(DNASequence(item["sequence"], name=item["name"]))
        
        to_fasta(export_objs, args.out) 
        print(f"Succès : {len(export_objs)} séquences écrites dans {args.out}")
    else:
        manual_index = 0
        for res in final_results:
            name = res.name if isinstance(res, DNASequence) else res["name"]
            sequence = res.sequence if isinstance(res, DNASequence) else res["sequence"]
            print(f">{name}_{manual_index}\n{sequence}") 
            manual_index += 1

def fast_a_food_splash():
    # ANSI Color Codes
    BUN    = "\033[38;5;214m" # Orange/Tan
    LETTUCE = "\033[32m"      # Green
    DNA    = "\033[36m"      # Cyan
    RESET  = "\033[0m"       # Reset to default
    BOLD   = "\033[1m"
    KETCHUP = "\033[31m"     # Red

    ascii_art = f"""
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
         {BOLD}"From PASTA to FASTA."{RESET}
    """
    print(ascii_art)


if __name__ == "__main__":

    main()
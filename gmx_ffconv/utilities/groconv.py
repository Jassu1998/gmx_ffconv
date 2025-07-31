from .utilities_groconv.coordinate_reorganise import read_gro_atoms, reorder_full_gro, reorder_full_gro_backconv, rewrite_gro_with_itp_data

def run_groconv(args):
    gro_data = read_gro_atoms(filename=args.coordfile)
    atom_lines = gro_data["atom_lines"]
    # Split molecule names and counts
    mol_names = args.name
    mol_counts = list(map(int, args.nmol))
    if len(mol_names) != len(mol_counts):
        raise ValueError("Number of molecule names and molecule counts must match.")
    molecules = list(zip(mol_names, mol_counts))
    reordered = reorder_full_gro(atom_lines, molecules, mapping_dir=args.mapping_dir)
    with open(args.output, 'w') as f:
        f.write(f"{gro_data['title']}\n")
        f.write(f"{gro_data['atom_count']}\n")
        f.writelines(reordered)
        f.write(f"{gro_data['box_line']}\n")
    print(f"Successfully generated {args.output}")
    if args.validate:
        gro_backconv = read_gro_atoms(filename=args.output)
        atom_lines_back = gro_backconv["atom_lines"]
        reordered_back = reorder_full_gro_backconv(atom_lines_back,molecules,mapping_dir=args.mapping_dir)
        with open(f"backconv_{args.coordfile}",'w') as f:
            f.write(f"{gro_backconv['title']}\n")
            f.write(f"{gro_backconv['atom_count']}\n")
            f.writelines(reordered_back)
            f.write(f"{gro_backconv['box_line']}\n")
        print(f"Successfully generated backconv_{args.coordfile}")


    if args.norename:
        print("ℹ️  Skipping renaming of atoms and residues (--norename).")
    else:
        gro_data = read_gro_atoms(f"{args.output}")
        fixed_atom_lines = rewrite_gro_with_itp_data(gro_data["atom_lines"], molecules, mapping_dir=args.mapping_dir)
        with open(args.output, 'w') as f:
            f.write(f"{gro_data['title']}\n")
            f.write(f"{gro_data['atom_count']}\n")
            f.writelines(fixed_atom_lines)
            f.write(f"{gro_data['box_line']}\n")
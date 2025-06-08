from .utilities_groconv.coordinate_reorganise import read_gro_atoms, reorder_full_gro

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
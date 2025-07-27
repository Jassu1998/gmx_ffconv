from .file_reader import read_atoms_section,read_bonds_section,read_settles_section
from .graph_functions import *


def run_ffmap_standard(args):
    atoms1 = read_atoms_section(args.itp1)
    bonds1 = read_bonds_section(args.itp1)
    if len(atoms1) != 1 and len(bonds1) == 0:
        bonds1 = read_settles_section(args.itp1)  # substitute bonds1 with settles
    atoms2 = read_atoms_section(args.itp2)
    bonds2 = read_bonds_section(args.itp2)
    if len(atoms2) != 1 and len(bonds2) == 0:
        bonds2 = read_settles_section(args.itp2)
    if args.duplicate:
        # Simple identity mapping: atom i → atom i
        mappings = [(i, i) for i in range(len(atoms1))]
        with open(f"mapping_{args.name}.csv", "w", newline='') as f:
            writer = csv.writer(f)
            # Write the header
            writer.writerow([args.itp1, args.itp2])
            # Write each mapping line
            for i, j in mappings:
                writer.writerow([i+1, j+1])
        with open(f"back_mapping_{args.name}.csv", "w", newline='') as f:
            writer = csv.writer(f)
            # Write the header
            writer.writerow([args.itp1, args.itp2])
            # Write each mapping line
            for i, j in mappings:
                writer.writerow([i+1, j+1])
        return
    if args.all_mappings and args.validate:
        raise ValueError("Cannot use --all-mappings and --validate together.")
    if args.all_mappings:
        mappings =match_graphs(atoms1, bonds1, atoms2, bonds2,all_mappings=True)
        mapping_writer(args, mappings=mappings)
    elif args.validate:
        mappings = match_graphs(atoms1, bonds1, atoms2, bonds2)
        mappings_back = match_graphs(atoms2, bonds2, atoms1, bonds1)
        mapping_writer(args,mappings=mappings,mapping_filename="mapping")
        mapping_writer(args,mappings=mappings_back,mapping_filename="back_mapping")
    else:
        mappings = match_graphs(atoms1, bonds1, atoms2, bonds2)  # all_mappings="false"
        mapping_writer(args,mappings=mappings)
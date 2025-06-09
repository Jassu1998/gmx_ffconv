from .file_reader import read_atoms_section,read_bonds_section
from .graph_functions import *

def run_ffmap_standard(args):
    atoms1 = read_atoms_section(args.itp1)
    if args.duplicate:
        # Simple identity mapping: atom i → atom i
        mappings = [(i, i) for i in range(len(atoms1))]
        for i in mappings:
            with open(f"mapping_{args.name}.csv", "w", newline='') as f:
                writer = csv.writer(f)
                writer.writerow([args.itp1, args.itp2])
                writer.writerow([i, i])
        print(f"Duplicated mapping written with {len(mappings)} identity pairs.")
        return
    bonds1 = read_bonds_section(args.itp1)
    atoms2 = read_atoms_section(args.itp2)
    bonds2 = read_bonds_section(args.itp2)
    if args.all_mappings:
        mappings =match_graphs(atoms1, bonds1, atoms2, bonds2,all_mappings=True)
    else:
        mappings = match_graphs(atoms1, bonds1, atoms2, bonds2)  # all_mappings="false"
    mapping_writer(args,mappings=mappings)
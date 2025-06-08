import os
import csv
def read_mapping(mapping_file):
    with open(mapping_file, newline='') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        return [(int(row[0]), int(row[1])) for row in reader]

def read_gro_atoms(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    title = lines[0].rstrip('\n')
    atom_count = int(lines[1].strip())
    atom_lines = lines[2:2 + atom_count]  # preserve \n
    box_line = lines[2 + atom_count].rstrip('\n')
    return {
        "title": title,
        "atom_count": atom_count,
        "atom_lines": atom_lines,
        "box_line": box_line
    }

def reorder_block_atoms(atom_lines, start_idx, mapping, n_mols):
    n_atoms_per_mol = len(mapping)
    expected_block_size = n_mols * n_atoms_per_mol
    if start_idx + expected_block_size > len(atom_lines):
        raise ValueError(f"Block exceeds available atom lines")
    reordered = [None] * expected_block_size
    for mol_index in range(n_mols):
        offset_old = start_idx + mol_index * n_atoms_per_mol
        offset_new = mol_index * n_atoms_per_mol
        for orig, new in mapping:
            orig_idx = offset_old + (orig - 1)
            new_idx = offset_new + (new - 1)
            reordered[new_idx] = atom_lines[orig_idx]
    if any(line is None for line in reordered):
        raise RuntimeError("Reordering failed. Mapping may be incomplete.")
    return reordered

def reorder_full_gro(atom_lines, molecules, mapping_dir="."):
    reordered_lines = []
    idx = 0
    for mol_name, mol_count in molecules:
        mapping_file = os.path.join(mapping_dir, f"mapping_{mol_name}.csv")
        if not os.path.isfile(mapping_file):
            raise FileNotFoundError(f"Mapping file not found: {mapping_file}")
        mapping = read_mapping(mapping_file)
        natoms_per_mol = len(mapping)
        for mol_index in range(mol_count):
            start = idx
            end = idx + natoms_per_mol
            if end > len(atom_lines):
                raise ValueError(f"Not enough lines for {mol_name} molecule {mol_index}")
            mol_lines = atom_lines[start:end]
            # Apply mapping to this individual molecule
            reordered_mol = [None] * natoms_per_mol
            for orig, new in mapping:
                reordered_mol[new - 1] = mol_lines[orig - 1]
            if any(x is None for x in reordered_mol):
                raise RuntimeError(f"Incomplete mapping for {mol_name} molecule {mol_index}")
            reordered_lines.extend(reordered_mol)
            idx += natoms_per_mol
    if idx != len(atom_lines):
        raise ValueError("Some atom lines were not processed. Mismatch in molecule counts?")
    return reordered_lines

def write_gro_file(filename, title, reordered_atoms, box_line):
    with open(filename, 'w') as f:
        f.write(f"{title}\n")
        f.write(f"{len(reordered_atoms)}\n")
        f.writelines(reordered_atoms)
        f.write(f"{box_line}\n")
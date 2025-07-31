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

def reorder_full_gro_backconv(atom_lines, molecules, mapping_dir="."):
    reordered_lines = []
    idx = 0
    for mol_name, mol_count in molecules:
        mapping_file = os.path.join(mapping_dir, f"back_mapping_{mol_name}.csv")
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


def get_itp_path_from_mapping(mapping_file):
    with open(mapping_file, newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        if len(header) < 2:
            raise ValueError(f"Mapping file {mapping_file} has fewer than two columns.")
        return header[1]  # Second column header: itp file path
def parse_itp_atoms(itp_path):
    atoms = []
    with open(itp_path) as f:
        lines = f.readlines()
    in_atoms = False
    for line in lines:
        stripped = line.strip()
        if not in_atoms:
            if stripped.lower().startswith('[ atoms'):
                in_atoms = True
        elif stripped == '' or stripped.startswith('['):
            break
        elif not stripped.startswith(';'):
            parts = stripped.split()
            if len(parts) >= 5:
                resname = parts[3]
                atomname = parts[4]
                atoms.append((resname, atomname))
    return atoms
def rewrite_gro_with_itp_data(atom_lines, molecules, mapping_dir="."):
    updated_lines = []
    idx = 0
    global_atom_id = 1
    global_res_id = 1

    def format_id(id_num):
        if id_num > 99999:
            return id_num % 100000
        else:
            return id_num

    mapping_cache = {}
    atom_defs_cache = {}

    for mol_name, mol_count in molecules:
        if mol_name not in mapping_cache:
            mapping_csv = os.path.join(mapping_dir, f"mapping_{mol_name}.csv")
            try:
                with open(mapping_csv) as f:
                    line = f.readline().strip()
                    parts = line.split(',')
                    if len(parts) < 2:
                        raise ValueError(f"Malformed line in {mapping_csv}: {line}")
                    itp_path_from_csv = parts[1].strip()
            except FileNotFoundError:
                raise FileNotFoundError(f"Mapping CSV not found: {mapping_csv}")
            except Exception as e:
                raise RuntimeError(f"Error reading {mapping_csv}: {e}")

            # Try direct path first
            if os.path.isfile(itp_path_from_csv):
                itp_file = itp_path_from_csv
            else:
                # Try relative to mapping_dir
                combined_path = os.path.join(mapping_dir, itp_path_from_csv)
                if os.path.isfile(combined_path):
                    itp_file = combined_path
                else:
                    # Prompt user
                    print(f"\n⚠️  ITP file for molecule '{mol_name}' not found.")
                    print(f"Tried: '{itp_path_from_csv}' and '{combined_path}'")
                    while True:
                        user_input = input(f"Please provide the full path to the ITP file for '{mol_name}': ").strip()
                        if os.path.isfile(user_input):
                            itp_file = user_input
                            break
                        else:
                            print(f"❌ File '{user_input}' does not exist. Please try again.")

            try:
                atom_defs = parse_itp_atoms(itp_file)
            except Exception as e:
                raise RuntimeError(f"Failed to parse ITP file '{itp_file}' for molecule '{mol_name}': {e}")

            mapping_cache[mol_name] = itp_file
            atom_defs_cache[mol_name] = atom_defs

        atom_defs = atom_defs_cache[mol_name]
        natoms_per_mol = len(atom_defs)
        for _ in range(mol_count):
            for i in range(natoms_per_mol):
                old_line = atom_lines[idx].rstrip('\n')
                resname, atomname = atom_defs[i]
                resname = resname[:5].ljust(5)
                atomname = atomname[:5].rjust(5)
                xyz = old_line[20:44]
                vel = old_line[44:]

                atom_id_fmt = format_id(global_atom_id)
                res_id_fmt = format_id(global_res_id)

                new_line = f"{res_id_fmt:5d}{resname}{atomname}{atom_id_fmt:5d}{xyz}{vel}"
                updated_lines.append(new_line + '\n')

                global_atom_id += 1
                idx += 1
            global_res_id += 1

    if idx != len(atom_lines):
        raise ValueError("Mismatch in total atom lines processed.")
    return updated_lines

def write_gro_file(filename, title, reordered_atoms, box_line):
    with open(filename, 'w') as f:
        f.write(f"{title}\n")
        f.write(f"{len(reordered_atoms)}\n")
        f.writelines(reordered_atoms)
        f.write(f"{box_line}\n")


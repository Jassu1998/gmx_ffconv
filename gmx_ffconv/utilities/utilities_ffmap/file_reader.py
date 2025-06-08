


def read_atoms_section(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    atoms = {}
    in_atoms = False
    for line in lines:
        stripped = line.strip()
        if not in_atoms:
            if stripped.lower().startswith('[ atoms'):
                in_atoms = True
        elif stripped == '':
            break
        elif not stripped.startswith(';'):
            parts = stripped.split()
            if len(parts) >= 5:
                atom_index = int(parts[0])
                atom_type = parts[4][0].upper() # First letter of atom type
                atoms[atom_index] = atom_type
    return atoms

def read_bonds_section(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    bonds = []
    in_bonds = False
    for line in lines:
        stripped = line.strip()
        if not in_bonds:
            if stripped.lower().startswith('[ bonds'):
                in_bonds = True
        elif stripped == '':
            break
        elif not stripped.startswith(';'):
            parts = stripped.split()
            if len(parts) >= 2:
                bonds.append((int(parts[0]), int(parts[1])))
    return bonds
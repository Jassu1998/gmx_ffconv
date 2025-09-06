from .file_reader import *
from .graph_functions import *
import csv
from collections import Counter
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher


def build_graph_naming(atoms, bonds):
    G = nx.Graph()
    for idx, attr in atoms.items():
        G.add_node(idx, props=[attr.get('atom_type'), attr.get('atom_name')])
    for i,j in bonds:
        G.add_edge(i,j)
    return G


def normalize_atoms(atoms1, atoms2, csv_atoms):
    # Step 1: Rename atoms1 according to CSV
    for attr in atoms1.values():
        if attr['atom_name'] in csv_atoms:
            attr['atom_name'] = csv_atoms[attr['atom_name']].strip()

    # Step 2: Build allowed names from renamed atoms1
    allowed_names = {
        attr['atom_name'].strip()
        for attr in atoms1.values()
        if attr['atom_name'] is not None
    }

    # Step 3: Clean atoms2 based on allowed names
    for attr in atoms2.values():
        name = attr['atom_name']
        if name is not None:
            name_clean = name.strip()
            if name_clean in allowed_names:
                attr['atom_name'] = name_clean
            else:
                attr['atom_name'] = None
    return allowed_names


def make_node_matcher(G1, G2, mode="relaxed"):
    def node_match(n1_attrs, n2_attrs, n1_id=None, n2_id=None):
        type1, name1 = n1_attrs.get('props', [None, None])
        type2, name2 = n2_attrs.get('props', [None, None])
        # exact name match
        if name1 and name2 and name1 == name2:
            return True
        # fallback: type must match
        if type1 != type2:
            return False
        # compare neighbors if node IDs are provided
        if n1_id is not None and n2_id is not None:
            neighbors1 = sorted([G1.nodes[n]['props'][0] for n in G1.neighbors(n1_id)])
            neighbors2 = sorted([G2.nodes[n]['props'][0] for n in G2.neighbors(n2_id)])
            return neighbors1 == neighbors2
        return True

    def strict_node_match(n1, n2):
        props1 = n1.get('props', [None, None])
        props2 = n2.get('props', [None, None])
        type1, name1 = props1
        type2, name2 = props2
        if name1 is not None and name2 is not None:
            return name1 == name2
        elif name1 is None and name2 is None:
            return type1 == type2
        else:
            return False

    return strict_node_match if mode == "strict" else node_match

def run_ffmap_naming(args):
    if args.validate:
        raise ValueError("The --validate option is not supported yet with consistent naming.")

    if args.duplicate:
        raise ValueError("The --duplicate option is not supported with consistent naming.")
    if args.all_mappings:
        raise ValueError("The option all_mappings is not supported yet with consistent naming.")

    atoms1 = read_atoms_section_atomname(args.itp1)
    bonds1 = read_bonds_section(args.itp1)
    if len(atoms1) != 1 and len(bonds1) == 0:
        bonds1 = read_settles_section(args.itp1)  # substitute bonds1 with settles
    atoms2 = read_atoms_section_atomname(args.itp2)
    bonds2 = read_bonds_section(args.itp2)
    if len(atoms2) != 1 and len(bonds2) == 0:
        bonds2 = read_settles_section(args.itp2)
    csv_atoms = {}
    with open(args.consistent_naming, newline='') as f:
        reader = csv.reader(f)
        next(reader)  # skip header if present
        for row in reader:
            atom_name = row[0].strip()  # first column
            new_name = row[1].strip()  # second column, optional
            csv_atoms[atom_name] = new_name
    # csv_atoms is your mapping from column1 -> column2
    allowed_names = set(csv_atoms.values())  # all names that are valid

    for idx, attr in atoms1.items(): # Rename atoms1 using CSV mapping
        if attr['atom_name'] in csv_atoms:
            attr['atom_name'] = csv_atoms[attr['atom_name']]

    for attr_dict in [atoms1, atoms2]: # After renaming, set any atom not in allowed_names to None
        for idx, attr in attr_dict.items():
            if attr['atom_name'] not in allowed_names:
                attr['atom_name'] = None

    # Sets collection of allowed names to those present int first topology
    allowed_names = set(attr['atom_name'].strip() for attr in atoms1.values() if attr['atom_name'] is not None)

    # Sets any atom names not present in the first topology to none
    for attr in atoms2.values():
        name = attr['atom_name']
        if name is not None:
            name_clean = name.strip()
            if name_clean in allowed_names:
                attr['atom_name'] = name_clean
            else:
                attr['atom_name'] = None
    # Get atom names, ignoring None
    atoms1_names = [attr['atom_name'] for attr in atoms1.values() if attr['atom_name'] is not None]
    atoms2_names = [attr['atom_name'] for attr in atoms2.values() if attr['atom_name'] is not None]
    # Count occurrences
    atoms1_count = Counter(atoms1_names)
    atoms2_count = Counter(atoms2_names)
    # Prints occurrences of each one
    print("Counts of each atom name in first topology:")
    for name, count in sorted(atoms1_count.items()):
        print(f"{name}: {count}")

    print("Counts of each atom name in the second topology:")
    for name, count in sorted(atoms2_count.items()):
        print(f"{name}: {count}")
    c1 = Counter(atoms1_names)
    c2 = Counter(atoms2_names)
    # Differences
    # Differences
    print("\nDifferences (if any):")
    all_names = set(c1.keys()) | set(c2.keys())  # union of all names
    for name in sorted(all_names):
        n1 = c1.get(name, 0)
        n2 = c2.get(name, 0)
        if n1 != n2:
            print(f"Atom name: {name}: first itp={n1}, second itp={n2}")



    G1 = build_graph_naming(atoms1, bonds1)
    G2 = build_graph_naming(atoms2, bonds2)

    if atoms1_count == atoms2_count:
        print("Counts match → using strict node matcher")
        matcher = GraphMatcher(G1, G2, node_match=make_node_matcher(G1, G2, mode="strict"))
    else:
        print("Counts differ → using relaxed node matcher. The performance will take a hit.")
        matcher = GraphMatcher(G1, G2, node_match=make_node_matcher(G1, G2, mode="relaxed"))

    if matcher.is_isomorphic():
        print("Graphs are isomorphic!")
        mappings = [matcher.mapping]
    else:
        # raise an error instead of silently computing something else
        raise ValueError("Graphs are NOT isomorphic – cannot produce a mapping")
    # Now 'mappings' exists in either case
    mapping_writer(args, mappings=mappings)
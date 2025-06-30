import networkx as nx
from networkx.algorithms import isomorphism
import csv
def build_graph(atoms, bonds):
    G = nx.Graph()
    for idx, atom_type in atoms.items():
        G.add_node(idx, atom_type=atom_type)
    for i, j in bonds:
        G.add_edge(i, j)
    return G

def match_graphs(atoms1, bonds1, atoms2, bonds2,all_mappings=False):
    G1 = build_graph(atoms1, bonds1)
    G2 = build_graph(atoms2, bonds2)
    nm = isomorphism.categorical_node_match('atom_type', None)
    matcher = isomorphism.GraphMatcher(G1, G2, node_match=nm)
    if all_mappings:
        mappings = []
        for mapping in matcher.isomorphisms_iter():
            mappings.append(mapping)
        if not mappings:
            raise ValueError("Graphs are not isomorphic.")
        return mappings
    else:
        try:
            mapping = next(matcher.isomorphisms_iter())
            return [mapping]  # return as a list for consistent interface
        except StopIteration:
            raise ValueError("Graphs are not isomorphic.")

def mapping_writer(args, mappings, mapping_filename="mapping"):
    # Decide header order
    if mapping_filename == "mapping":
        header = [args.itp1, args.itp2]
    elif mapping_filename == "back_mapping":
        header = [args.itp2, args.itp1]
    else:
        header = [args.itp1, args.itp2]  # default

    if len(mappings) == 1:
        mapping = mappings[0]
        matched_indices = sorted(mapping.items(), key=lambda x: x[1])
        filename = f"{mapping_filename}_{args.name}.csv"
        with open(filename, "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(matched_indices)
    else:
        for i, mapping in enumerate(mappings, start=1):
            matched_indices = sorted(mapping.items(), key=lambda x: x[1])
            filename = f"{mapping_filename}_{args.name}_{i}.csv"
            with open(filename, "w", newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerows(matched_indices)


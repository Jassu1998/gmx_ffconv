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

def mapping_writer(args, mappings):
    if len(mappings) == 1:
        # single file
        mapping = mappings[0]
        matched_indices = sorted(mapping.items(), key=lambda x: x[1])
        back_matched = sorted(mapping.items(), key=lambda x: x[0])
        back_match_sorted = [(b, a) for (a, b) in back_matched]
        with open(f"mapping_{args.name}.csv", "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow([args.itp1, args.itp2])
            writer.writerows(matched_indices)
        with open(f"back_mapping_{args.name}.csv", "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow([args.itp2, args.itp1])
            writer.writerows(back_match_sorted)

    else:
        # multiple files
        for i, mapping in enumerate(mappings, start=1):
            matched_indices = sorted(mapping.items(), key=lambda x: x[1])
            back_matched = sorted(mapping.items(), key=lambda x: x[0])
            back_match_sorted = [(b, a) for (a, b) in back_matched]
            with open(f"mapping_{args.name}_{i}.csv", "w", newline='') as f:
                writer = csv.writer(f)
                writer.writerow([args.itp1, args.itp2])
                writer.writerows(matched_indices)
            with open(f"back_mapping_{args.name}_{i}.csv", "w", newline='') as f:
                writer = csv.writer(f)
                writer.writerow([args.itp2, args.itp1])
                writer.writerows(back_match_sorted)

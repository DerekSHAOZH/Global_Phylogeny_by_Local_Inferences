from ete3 import Tree
import sys

def parse_newick_string(filename_our_method, filename_benchmark):
    newick_our_method = ""
    newick_raxml_benchmark = ""
    with open(filename_our_method, 'r') as file_our_method:
        newick_our_method = file_our_method.readlines()[-1]
    with open(filename_benchmark, 'r') as file_benchmark:
        newick_raxml_benchmark = file_benchmark.read()
    return newick_our_method, newick_raxml_benchmark

def calculate_rf_distance(newick_our_method, newick_benchmark):
    tree_our_method = Tree(newick_our_method)
    tree_benchmark = Tree(newick_benchmark)
    (rf, rf_max, names, edges_our_method, edges_benchmark, discarded_edges_our_method, discarded_edges_benchmark) = tree_our_method.robinson_foulds(tree_benchmark, unrooted_trees=True)
    edges_union = edges_our_method.union(edges_benchmark)
    edges_intersection = edges_our_method.intersection(edges_benchmark)

    # RF(T1, T2) = |S1⋃S2| − |S1⋂S2|
    print(len(edges_our_method))
    RF_dist = edges_union - edges_intersection
    print(f'Robinson-Foulds distance = |S1⋃S2| − |S1⋂S2| = {len(edges_union)} - {len(edges_intersection)} = {len(edges_union) - len(edges_intersection)} over a total of {rf_max}')
    return RF_dist

if __name__ == "__main__":
    newick_our_method, newick_raxml_benchmark = parse_newick_string(sys.argv[1], sys.argv[2])
    calculate_rf_distance(newick_our_method, newick_raxml_benchmark)
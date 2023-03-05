# import config
from config import *
from mst import *
from tree import *
#I. Selecting Vs from the MST such that it induces a subtree
def find_subtree(T, threshold): # e.g. threshold = 10
  # initialize vertices_to_visit to all of the leaves
  vertices_to_visit = T.leaves
  k = 1
  while len(vertices_to_visit) > 0:
    #randomize the starting vertex
    if k == 1:
      v = vertices_to_visit.pop(random.randrange(len(vertices_to_visit)))
      # print(v.name)
    else:
      v = vertices_to_visit.pop(0)
  
    if v.is_leaf:
        v.subtree.append(v)
        if v.neighbors[0] not in vertices_to_visit:
          #randomly shuffle the neighbors
          vertices_to_visit.extend(random.sample(v.neighbors, len(v.neighbors)))
    else:
    #if v is an internal vertex, 
      #check the number of neighbors with empty subtrees
      num_neighbor_with_empty_subtree = sum([len(n.subtree) == 0 for n in v.neighbors])
      # print(num_neighbor_with_empty_subtree)

      #if only one of its neighbors has an empty(unassigned) subtree, assign the subtree at this node
      if num_neighbor_with_empty_subtree == 1:
        # print(v.name)
        v.subtree.append(v)
        #randomly shuffle the neighbors
        for n in random.sample(v.neighbors, len(v.neighbors)):
          if len(n.subtree) > 0:
            for subtree_vertex in n.subtree:
              if subtree_vertex not in v.subtree:
                v.subtree.append(subtree_vertex)
          elif n not in vertices_to_visit:
            vertices_to_visit.append(n)

      #if the number of neighbors with empty(unassigned) subtrees > 1, 
      #then queue this vertex for re-visiting after updating its neighbors
      else:
        #randomly shuffle the neighbors
        for n in random.sample(v.neighbors, len(v.neighbors)):
          if len(n.subtree) == 0 and n not in vertices_to_visit:
            vertices_to_visit.append(n)
        vertices_to_visit.append(v)
 
    if len(v.subtree) >= threshold: 
      return v.subtree

    k += 1

#II. Selecting Ve by performing BFS from the root of the induced subtree from Vs 
#input: MST tree, Vs, threshold for Ve which is consistent with Vs
#output: Ve consisting of neighboring vertices of Vs
def extra_vertices_with_bfs(MST, Vs, threshold): #function for BFS
  Ve = [] 
  visited = []  #List for visited nodes.
  queue = []     #Initialize a queue

  visited.append(Vs[0])
  queue.append(Vs[0])

  while queue:          # Creating loop to visit each node
    m = queue.pop(0) 
    
    for neighbor in random.sample(MST.Get_vertex(m.name).neighbors, len(MST.Get_vertex(m.name).neighbors)):
      if neighbor not in visited:
        visited.append(neighbor)
        queue.append(neighbor)
        if neighbor not in Ve and neighbor not in Vs: #vertices in Vs should not be in Ve
          Ve.append(neighbor)
      if len(Ve) >= threshold:
        return Ve

# #III. Finding all internal edges for NNI
# def find_internal_edge(T):
#   internal_vertex_names = [v.name for v in T.vertex_map.values() if v.is_leaf == False]
#   internal_edge_possibilities = list(itertools.combinations(internal_vertex_names,2))
#   internal_edges = [p for p in internal_edge_possibilities if T.Get_vertex(p[0]) in T.Get_vertex(p[1]).neighbors]
#   return internal_edges

# #IV. Searching for all possible tree topologies with NNI
# #input: unrooted tree
# def search_tree(tree, root_branch_node1, root_branch_node2, nni_node1, nni_node2):
#   # 1. COMPUTE A TREE
#   # root the tree at an arbitrary branch
#   tree.Root_tree_along_branch(root_branch_node1, root_branch_node2)

#   # 2. COMPUTE THE PARSIMONY SCORE (of the rooted tree)
#   tree.Run_Sankoff_for_all_sites()
#   total_parsimony_score = tree.total_parsimony_score
#   # print(f"Total parsimony score is {tree.total_parsimony_score}")

#   # 3. MODIFY TREE TOPOLOGY
#   tree.Suppress_the_root()
#   T1, T2 = tree.NNI(nni_node1, nni_node2)

#   T1.Root_tree_along_branch(nni_node1, nni_node2)
#   T2.Root_tree_along_branch(nni_node1, nni_node2)
#   T1.Run_Sankoff_for_all_sites()
#   T2.Run_Sankoff_for_all_sites()
#   T1_parsimony_score = T1.total_parsimony_score
#   T2_parsimony_score = T2.total_parsimony_score
#   # print(f"T1 parsimony score: {T1_parsimony_score}")
#   # print(f"T2 parsimony score: {T2_parsimony_score}")

#   if T1_parsimony_score < total_parsimony_score and T1_parsimony_score < T2_parsimony_score:
#     # print("changing to T1")
#     T1.Suppress_the_root()
#     return T1
#   elif T2_parsimony_score < total_parsimony_score and T2_parsimony_score < T1_parsimony_score:
#     # print("changing to T2")
#     T2.Suppress_the_root()
#     return T2
#   else:
#     # print("not changing tree topology")
#     return tree   


#V. Converting sequence dictionary to alignment object
def convert_seq_to_aln(V, sequences, mapping):

  sequences_sub = {mapping[k]: sequences[k] for k in [v.name for v in V]}
  count = len(sequences_sub)
  length = max(map(len, sequences_sub.values()))
  msa = f" {count} {length}\n"
  msa += '\n'.join(f"{prot_id:<10} {sequence}" for prot_id, sequence in sequences_sub.items())
  # print(msa) 
  aln = AlignIO.read(StringIO(msa), 'phylip')
  return aln

#VI. Converting MP tree object to newick format
def newick_MP_local_tree(tree, mapping):
  Phylo.write(tree, "tmp.tre", 'newick')

  with open("tmp.tre", 'r') as f:
    output = f.read()
  tmp = output[0:-1]
  t = ete3.Tree(tmp, format = 1)
  newick_str = t.write(format = 9)
  for v_name, index in mapping.items():
    newick_str = newick_str.replace(index, v_name)
    # print(newick_str)
  return newick_str

#VII. Finding the local tree with maximum parsimony score
def infer_MP_local_tree(V, sequences, iteration_counter, cost_matrix):
  #The sequence ID for Phylo Alignment object can be up to 10 characters long
  #So, create a dictionary mapping V name to some shorter indices 
  mapping = {}
  for i in range(len(V)):
    index = 'No' +str(i) +'de'
    mapping[V[i].name] = index
  # print(mapping)
  aln = convert_seq_to_aln(V, sequences, mapping)

  # final_MP_newick = ''
  # final_lowest_parsimony_score = float('inf')

  # #starting with 3 random initial tree topologies to avoid local optimum
  # k = 1
  # while k < 4:
  #random starting tree
  t = ete3.Tree()
  t.populate(len(mapping), [v for v in mapping.values()])
  newick_str = t.write(format = 3)
  # print(newick_str)
  starting_tree = Phylo.read(StringIO(newick_str), "newick")
  # print(starting_tree)
  scorer = ParsimonyScorer()
  searcher = NNITreeSearcher(scorer)
  constructor = ParsimonyTreeConstructor(searcher, starting_tree)
  MP_tree = constructor.build_tree(aln)

  # cur_lowest_parsimony_score = scorer.get_score(MP_tree, aln)
  lowest_parsimony_score = scorer.get_score(MP_tree, aln)
  MP_newick = newick_MP_local_tree(MP_tree, mapping)
  # print(MP_newick)
    # if cur_lowest_parsimony_score < final_lowest_parsimony_score:
    #   final_lowest_parsimony_score = cur_lowest_parsimony_score
    #   final_MP_newick = MP_newick

  # print(f'Best parsimony score: {final_lowest_parsimony_score}')
  print(f'Best parsimony score: {lowest_parsimony_score}')
  T = Tree('local')
  # T.Read_newick_string_without_branch_lengths(final_MP_newick, iteration_counter)
  T.Read_newick_string_without_branch_lengths(MP_newick, iteration_counter)
  for v in V:
    # print(v.name)
    T.Get_vertex(v.name).sequence = sequences[v.name]

  T.sequence_length = len(sequences[V[0].name])
  
  T.Add_mut_cost_matrix(cost_matrix)
  T.Store_pos_list_for_unique_site_patterns()
  # print(f"Total number of unique site patterns is {len(T.first_pos_of_unique_site_pattern_to_pos_list)}")
  T.Get_total_number_of_informative_site_patterns()
  T.Run_Sankoff_for_all_sites()
  # print(f"Total parsimony score is {T.total_parsimony_score}")
  T.Suppress_the_root()
  # T.total_parsimony_score = lowest_parsimony_score

  return T

# #VII. Finding the local tree with maximum parsimony score
# def infer_MP_local_tree(V, sequences, iteration_counter, cost_matrix):

#   final_lowest_parsimony_score = float('inf')
#   k = 1

#   #starting with 3 random initial tree topologies to avoid local optimum
#   while k < 4:
#     print(f'##### Infer MP local tree: interation {k} begins #####')
#     t = ete3.Tree()
#     t.populate(len(V), [v.name for v in V])
#     newick_str = t.write(format = 9)
    
    
#     T = Tree('local')
#     T.Read_newick_string_without_branch_lengths(newick_str, iteration_counter)
#     for v in V:
#       T.Get_vertex(v.name).sequence = sequences[v.name]
#       #print(T.Get_vertex(v.name).sequence)
#     T.sequence_length = len(sequences[V[0].name])
#     # for l in T.leaves:
#     #   print(l.sequence)

#     # for v in T.vertex_map.values():
#     #   print(f'vertex {v.name}: out-degree {v.out_degree}, neighbors {[n.name for n in v.neighbors]}, parent {v.parent.name}, children {[c.name for c in v.children]}')

#     T.Add_mut_cost_matrix(cost_matrix)
#     T.Store_pos_list_for_unique_site_patterns()
#     print(f"Total number of unique site patterns is {len(T.first_pos_of_unique_site_pattern_to_pos_list)}")
#     T.Get_total_number_of_informative_site_patterns()
#     T.Run_Sankoff_for_all_sites()
#     total_parsimony_score = T.total_parsimony_score
#     print(f"Total parsimony score is {T.total_parsimony_score}")
#     T.Suppress_the_root()

#     if k == 1:
#       T_final = T

#     T_lowest = T
#     cur_lowest_parsimony_score = float('inf')
#     continue_search = True

#     while(continue_search):
#       # print("one more while iteration")
#       T_new = T_lowest
#       continue_search = False
#       internal_edges = find_internal_edge(T_lowest)
      

#       for b in internal_edges:
#         # print(f'NNI branch: ({b[0]}, {b[1]})')
#         tree = search_tree(T_new, b[0], b[1], b[0], b[1])
#         if tree.total_parsimony_score < cur_lowest_parsimony_score:
#           # print("Update continue_search to TRUE")
#           cur_lowest_parsimony_score = tree.total_parsimony_score
#           T_lowest = tree
#           continue_search = True
#     print(f'For iteration {k}, best parsimony score: {cur_lowest_parsimony_score}')

#     if cur_lowest_parsimony_score < final_lowest_parsimony_score:
#       T_final = T_lowest
#       final_lowest_parsimony_score = cur_lowest_parsimony_score

#     k += 1

#   return T_final

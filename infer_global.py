# import config
from config import *
from infer_local import *
from mst import *
from tree import *
#I. Finding the subtree induced by only Vs in the local tree
#input: local unrooted phylogeny tree, Vs
#output: dictionary with a tuple of two child nodes from Vs as key and the corresponding hidden vertex as value 
#e.g.,    A
#        / \
#       B   C
def informative_subtree(tree, Vs):
  forest = {}
  for v_tuple in itertools.combinations(Vs, 2):
    if tree.Get_vertex(v_tuple[0].name).neighbors[0] is tree.Get_vertex(v_tuple[1].name).neighbors[0]:
      forest[v_tuple] = tree.Get_vertex(v_tuple[0].name).neighbors[0]
  
  return forest

#II. Deleting paired Vs vertices and Adding the new parent vertices to V
def update_subtree_vertices(forest, V):
  V_new = V.copy()
  
  for children, parent_new in forest.items():
    if children[0] in V:
      V_new.remove(children[0])
      # print(f'removed {children[0].name} from V_new')
    if children[1] in V:
      V_new.remove(children[1])
      # print(f'removed {children[1].name} from V_new')
    V_new.append(parent_new)
    # print(f'added {parent_new.name} to V_new')
  return V_new

#III. Deleting paired Vs vertices and Adding the new parent vertices to the sequence dictionary
def update_sequences(forest, seq_dict):
  sequences_new = seq_dict.copy()
  for children, parent_new in forest.items():
    if children[0].name in seq_dict.keys():
      del sequences_new[children[0].name]
      # print(f'removed {children[0].name} from sequences_new')
    if children[1].name in seq_dict.keys():
      del sequences_new[children[1].name]
      # print(f'removed {children[1].name} from sequences_new')
    sequences_new[parent_new.name] = parent_new.sequence
    # print(f'added {parent_new.name} to sequences_new')
  return sequences_new


#IV. Updating MST
def update_mst(MST, forest, V_old, sequences_sub, dist_mat_sub):
  MST_copy = MST.copy()
  #delete all within-subgraph edge from MST
  M_sub = nx.subgraph(MST, [v.name for v in V_old])
  for edge in list(M_sub.edges):
    MST_copy.remove_edge(edge[0], edge[1])
  #delete paired $V_s$ vertices and add new parent vertices to MST
  for children, parent_new in forest.items():
      if children[0].name in MST_copy.nodes:
        MST_copy.remove_node(children[0].name)
        # print(f'removed {children[0].name} from MST')
      if children[1].name in MST_copy.nodes:
        MST_copy.remove_node(children[1].name)
        # print(f'removed {children[1].name} from MST')
      MST_copy.add_node(parent_new.name)
      # print(f'added {parent_new.name} to MST')
  #create the sub-MST of V
  MST_sub = create_mst_from_dist_mat(sequences_sub, dist_mat_sub)
  #add edges in sub-MST to MST
  for edge in list(MST_sub.edges):
    MST_copy.add_edge(edge[0], edge[1], weight = MST_sub.get_edge_data(edge[0], edge[1])['weight'])

  return MST_copy  

#V. Inferring the global phylogeny tree
def infer_global_tree(T_global, MST, sequences, cost_matrix, threshold):
  #k is an int as the iteration counter to avoid duplicated naming for hidden vertex
  k = 1
  #is_last is a boolean variable indicating whether this is the last iteration
  is_last = False
  #num_sequences is an int storing the initial number of sequences (taxons)
  num_sequences = len(MST.nodes)
  #plot_1 is a boolean variable indicating whether to plot MST when number of nodes drops to 2/3 of the initial number for the first time
  plot_1 = True
  #plot_2 is a boolean variable indicating whether to plot MST when number of nodes drops to 1/3 of the initial number for the first time
  plot_2 = True

  while len(T_global.edge_list_map) < len(T_global.vertex_map) - 1:
    print(f'#################### Infer global tree: iteration {k} begins ####################')

    #Converting networkx object to our own Tree class
    M = Tree('MST')
    for edge_tuple in MST.edges():
      M.Add_edge(edge_tuple[0], edge_tuple[1], MST.get_edge_data(edge_tuple[0], edge_tuple[1])['weight'])
    M.Set_MST_leaf_flags()
    M.Set_leaves()
    for vertex_name in M.vertex_map.keys():
      M.Get_vertex(vertex_name).subtree = []

    #if number of vertices in MST is smaller than threshold,
    if len(M.vertex_map) <= threshold:
      #then all vertices will be used to constitute the last local tree
      V = [i for i in M.vertex_map.values()]
      is_last = True
    
    else:
      #Selecting Vs from the MST such that it induces a subtree
      print('Finding Vs ...')
      Vs = find_subtree(M, threshold)

      #if |Vm| - |Vs| > threshold,
      if(len(M.vertex_map) - len(Vs) > threshold):
        #Selecting Ve by performing BFS from the root of the induced subtree from Vs
        print('Finding Ve ...')
        Ve = extra_vertices_with_bfs(M, Vs, threshold)
          
        V = Vs + Ve
      else:
        V = [i for i in M.vertex_map.values()]
        is_last = True
    
    print('Inferring MP local tree ...')
    T_lowest = infer_MP_local_tree(V, sequences, k, cost_matrix)

    #if this is the last iteration, 
    if is_last:
      #add all edges in the local tree to the global tree
      vertices_visited = []
      for v_name, v in T_lowest.vertex_map.items():
        vertices_visited.append(v_name)
        for neighbor in v.neighbors:
          if(neighbor.name not in vertices_visited):
            T_global.Add_edge(v_name, neighbor.name, 0.1)
            print(f'Updating global tree: Added edge from {v_name} to {neighbor.name}, distance = 0.1')
        if T_global.Get_vertex(v_name).sequence != '':
          T_global.Get_vertex(v_name).sequence = v.sequence
      
    else: 
      #Finding the subtree induced by only Vs in the local tree
      forest = informative_subtree(T_lowest, Vs)
#       print(f'Forest: {len(forest)}')

      #if forest is not empty, 
      if (len(forest) > 0):

        #Updating the global phylogeny tree
        print('Updating global tree ...')
        for children, parent_new in forest.items():
          T_global.Add_edge(children[0].name, parent_new.name, 0.1)
          T_global.Add_edge(children[1].name, parent_new.name, 0.1)
          print(f'Updating global tree: Added edge from {children[0].name} to {parent_new.name}, distance = 0.1')
          print(f'Updating global tree: Added edge from {children[1].name} to {parent_new.name}, distance = 0.1')
          T_global.Get_vertex(parent_new.name).sequence = parent_new.sequence


        #Deleting paired Vs vertices and Adding the new parent vertices to V
        V_new = update_subtree_vertices(forest, V)

        #Deleting paired Vs vertices and Adding the new parent vertices to the sequence dictionary
        sequences = update_sequences(forest, sequences)
        
        #Updating the subgraph distance matrix of V_new
        sequences_sub = {k: sequences[k] for k in [v.name for v in V_new]}
        dist_mat_sub = create_distance_matrix(sequences_sub)
        
        #Updating the edges in the subgraph consisting of V_new in MST
        print('Updating MST ...')
        MST = update_mst(MST, forest, V, sequences_sub, dist_mat_sub)

      print(f'Current MST has {len(MST.nodes)} nodes')
      if plot_1 and len(MST.nodes) <= (2/3)*num_sequences:
        fig, ax = plt.subplots(figsize=(30, 30))
        nx.draw(MST, with_labels=True, font_weight='bold', node_size = 500)
        plot1_fn = fn + f'.MST_{len(MST.nodes)}.png'
        fig.savefig(plot1_fn)
        plot_1 = False
      if plot_2 and len(MST.nodes) <= (1/3)*num_sequences:
        fig, ax = plt.subplots(figsize=(30, 30))
        nx.draw(MST, with_labels=True, font_weight='bold', node_size = 500)
        plot2_fn = fn + f'.MST_{len(MST.nodes)}.png'
        fig.savefig(plot2_fn)
        plot_2 = False

      print(f'Current global tree has {len(T_global.vertex_map)} vertices and {len(T_global.edge_list_map)} edges')
    
    

    df.loc[k] = [k, T_lowest.total_parsimony_score, len(MST.nodes), len(T_global.vertex_map), len(T_global.edge_list_map), len(V)]

    k += 1
  
  return T_global

def newick_global_tree(tree):
  last_iter = '_' + str(max(df['k']))
  for edge_tuple in tree.edge_list_map.keys():
    if last_iter in edge_tuple[0].name and last_iter in edge_tuple[1].name:
      branch = edge_tuple
      break

  #Add 1 to out_degree of all hidden vertices
  for v_name, v in tree.vertex_map.items():
    if 'EPI' not in v_name:
      v.out_degree = 1

  tree.Root_tree_along_branch(branch[0].name, branch[1].name)
  newick_str = tree.Compute_newick_format()
  return newick_str

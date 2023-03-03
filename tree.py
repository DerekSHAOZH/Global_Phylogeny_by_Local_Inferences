#I. Importing the Tree Class for phylogeny inference
class Vertex:
  def __init__(self,name): 
    self.name = name
    self.in_degree = 0
    self.out_degree = 0
    self.parent = self
    self.children = []
    self.neighbors = []
    self.newick_label = ""
    self.sequence = ""
    self.min_cost_subtree = []
    self.is_leaf = False
    self.is_root = False
    self.subtree = []

class Tree:
  def __init__(self,name):
    self.name = name
    self.vertex_map = {}
    self.pre_order_list = []
    self.post_order_list = []
    self.edge_list_map = {}
    self.root = -1
    self.mut_cost_matrix = np.array([[0, 1, 1, 1], 
                                 [1, 0, 1, 1],
                                 [1, 1, 0, 1],
                                 [1, 1, 1, 0]]) 
    self.total_parsimony_score = 0
    self.DNA_ind_map = {"A":0,"T":1,"G":2,"C":3}
    self.DNA_list = "ATGC"
    self.max_penalty = float("inf")    
    self.first_pos_of_unique_site_pattern_to_pos_list = {}
    self.leaves = []
    self.sequence_length = 0
  def Add_mut_cost_matrix(self,W):
    self.mut_cost_matrix = W  
  def Add_vertex(self,name):    
    v = Vertex(name)
    self.vertex_map[name] = v  
  def Contains_vertex(self,name):
    return (name in self.vertex_map.keys())
  def Get_vertex(self,name):
    if name in self.vertex_map.keys():
      return (self.vertex_map[name])

  def Add_edge(self, end1_name, end2_name, distance):   #FOR PROJECT
    if end1_name not in self.vertex_map.keys():
      self.Add_vertex(end1_name)
    if end2_name not in self.vertex_map.keys():
      self.Add_vertex(end2_name)
    node1 = self.Get_vertex(end1_name)
    node2 = self.Get_vertex(end2_name)   
    node1.neighbors.append(node2)
    node2.neighbors.append(node1)
    self.edge_list_map[(node1,node2)] = distance

  def Add_directed_edge(self, parent_name, child_name, distance):   
    if parent_name not in self.vertex_map.keys():
      self.Add_vertex(parent_name)
    if child_name not in self.vertex_map.keys():
      self.Add_vertex(child_name)
    p = self.Get_vertex(parent_name)
    c = self.Get_vertex(child_name)    
    p.out_degree += 1
    c.in_degree += 1
    c.parent = p
    p.children.append(c)
    self.edge_list_map[(p,c)] = distance
  def Get_edge_length(self, parent, child):
    return (self.edge_list_map[(parent, child)])
  def Set_root(self):
    for vertex in self.vertex_map.values():
      if vertex.in_degree == 0:
        self.root = vertex
        self.root.is_root = True
      else:
        vertex.is_root = False
  def Get_root(self):
    if self.root == -1:
      self.Set_root()
    return (self.root)

  def Set_MST_leaf_flags(self):   #FOR PROJECT
    for v in self.vertex_map.values():
      if len(v.neighbors) == 1:
        v.is_leaf = True
      else:
        v.is_leaf = False
        
  def Set_leaf_flags(self):
    for v in self.vertex_map.values():
      if v.out_degree == 0:
        v.is_leaf = True
      else:
        v.is_leaf = False
  def Set_leaves(self):
    self.leaves = []
    for v in self.vertex_map.values():      
      if v.is_leaf:
        self.leaves.append(v)    
  def Set_pre_order_and_post_order(self):
    self.Set_root()
    self.Set_leaf_flags()
    self.Set_leaves()
    self.pre_order_list = [self.root]
    self.post_order_list = [self.root]
    vertices_to_visit = [self.root]
    while len(vertices_to_visit) > 0:
      v = vertices_to_visit.pop()
      vertices_to_visit += v.children      
      self.post_order_list = v.children + self.post_order_list
      self.pre_order_list = self.pre_order_list + v.children
  def Compute_newick_format(self):
    if len(self.post_order_list) != len(self.vertex_map):
      self.Set_pre_order_and_post_order()
    for v in self.post_order_list:
      if v.out_degree == 0:
        v.newick_label = v.name
      else:        
        c_l = v.children[0]
        c_r = v.children[1]        
        len_l = self.Get_edge_length(v,c_l)
        len_r = self.Get_edge_length(v,c_r)
        v.newick_label = f'({c_l.newick_label}:{len_l},{c_r.newick_label}:{len_r})'
    self.root.newick_label += ";"
    return(self.root.newick_label)   

  ################################################
  # Add an iteration counter to avoind duplicated names for hidden vertices
  #e.g., hidden_vertex_name = "h" + str(hidden_vertex_ind) + '_' + str(k)
  ################################################
  def Read_newick_string(self, newick_string, iteration_num):    
    rx = r'\([^()]+\)'
    hidden_vertex_ind = 1
    while "," in newick_string:                  
      # search for the parenthesis
      m = re.search(rx,newick_string)
      # returns a tuple containing all the subgroups of the match "()"
      string_match = m.group()            
      # remove ( and )
      siblings_string = string_match[1:-1]      
      c_left_name_and_length, c_right_name_and_length = siblings_string.split(",")
      c_left_name, c_left_length = c_left_name_and_length.split(":")
      c_right_name, c_right_length = c_right_name_and_length.split(":")
      if not self.Contains_vertex(c_left_name):
          self.Add_vertex(c_left_name)
      if not self.Contains_vertex(c_right_name):
          self.Add_vertex(c_right_name)
      hidden_vertex_name = "h" + str(hidden_vertex_ind) + '_' + str(iteration_num)
      self.Add_vertex(hidden_vertex_name)            
      self.Add_directed_edge(hidden_vertex_name, c_left_name, float(c_left_length))
      self.Add_directed_edge(hidden_vertex_name, c_right_name, float(c_right_length))
      newick_string = newick_string.replace(string_match,hidden_vertex_name)
      hidden_vertex_ind += 1 

  def Read_newick_string_without_branch_lengths(self,newick_string, iteration_num):
    rx = r'\([^()]+\)'
    hidden_vertex_ind = 1
    while "," in newick_string:                  
      # search for the parenthesis
      m = re.search(rx,newick_string)
      # returns a tuple containing all the subgroups of the match "()"
      string_match = m.group()            
      # remove ( and )
      siblings_string = string_match[1:-1]      
      c_left_name, c_right_name = siblings_string.split(",")
      # c_left_name, c_left_length = c_left_name_and_length.split(":")
      # c_right_name, c_right_length = c_right_name_and_length.split(":")
      if not self.Contains_vertex(c_left_name):
          self.Add_vertex(c_left_name)
      if not self.Contains_vertex(c_right_name):
          self.Add_vertex(c_right_name)
      hidden_vertex_name = "h" + str(hidden_vertex_ind) + '_' + str(iteration_num)
      self.Add_vertex(hidden_vertex_name)            
      self.Add_directed_edge(hidden_vertex_name, c_left_name, 0.001)
      self.Add_directed_edge(hidden_vertex_name, c_right_name, 0.001)      
      newick_string = newick_string.replace(string_match,hidden_vertex_name)
      hidden_vertex_ind += 1 

  # def Read_newick_string(self, newick_string):    
  #   rx = r'\([^()]+\)'
  #   hidden_vertex_ind = 1
  #   while "," in newick_string:                  
  #     # search for the parenthesis
  #     m = re.search(rx,newick_string)
  #     # returns a tuple containing all the subgroups of the match "()"
  #     string_match = m.group()            
  #     # remove ( and )
  #     siblings_string = string_match[1:-1]      
  #     c_left_name_and_length, c_right_name_and_length = siblings_string.split(",")
  #     c_left_name, c_left_length = c_left_name_and_length.split(":")
  #     c_right_name, c_right_length = c_right_name_and_length.split(":")
  #     if not self.Contains_vertex(c_left_name):
  #         self.Add_vertex(c_left_name)
  #     if not self.Contains_vertex(c_right_name):
  #         self.Add_vertex(c_right_name)
  #     hidden_vertex_name = "h" + str(hidden_vertex_ind)
  #     self.Add_vertex(hidden_vertex_name)            
  #     self.Add_directed_edge(hidden_vertex_name, c_left_name, float(c_left_length))
  #     self.Add_directed_edge(hidden_vertex_name, c_right_name, float(c_right_length))
  #     newick_string = newick_string.replace(string_match,hidden_vertex_name)
  #     hidden_vertex_ind += 1 


  # def Read_newick_string_without_branch_lengths(self,newick_string):
  #   rx = r'\([^()]+\)'
  #   hidden_vertex_ind = 1
  #   while "," in newick_string:                  
  #     # search for the parenthesis
  #     m = re.search(rx,newick_string)
  #     # returns a tuple containing all the subgroups of the match "()"
  #     string_match = m.group()            
  #     # remove ( and )
  #     siblings_string = string_match[1:-1]      
  #     c_left_name, c_right_name = siblings_string.split(",")
  #     # c_left_name, c_left_length = c_left_name_and_length.split(":")
  #     # c_right_name, c_right_length = c_right_name_and_length.split(":")
  #     if not self.Contains_vertex(c_left_name):
  #         self.Add_vertex(c_left_name)
  #     if not self.Contains_vertex(c_right_name):
  #         self.Add_vertex(c_right_name)
  #     hidden_vertex_name = "h" + str(hidden_vertex_ind)
  #     self.Add_vertex(hidden_vertex_name)            
  #     self.Add_directed_edge(hidden_vertex_name, c_left_name, 0.001)
  #     self.Add_directed_edge(hidden_vertex_name, c_right_name, 0.001)      
  #     newick_string = newick_string.replace(string_match,hidden_vertex_name)
  #     hidden_vertex_ind += 1 


  def Read_fasta_file(self,fasta_file_name):
    fastaFile = open(fasta_file_name,'r')
    seq = ''
    name = ''
    sequenceAlignment={}
    for line in fastaFile:
        if line.startswith('>'):
            if seq != '':
                seq = seq.upper()
                if not self.Contains_vertex(name):
                  self.Add_vertex(name)
                v = self.Get_vertex(name)
                v.sequence = seq       
                seq = ''
            name = line.strip().split('>')[1]
        else:
            seq += line.strip()
    if not self.Contains_vertex(name):
      self.Add_vertex(name)
    v = self.Get_vertex(name)
    v.sequence = seq    
    self.sequence_length = len(seq)
    # sequenceAlignment[name] = seq
    fastaFile.close()
  def Set_null_sequences_in_ancestors(self):
    for v in self.vertex_map.values():
      if not v.is_leaf:
        v.sequence = "N"*self.sequence_length
  def Get_site_pattern(self,site):
    site_pattern = "".join([l.sequence[site] for l in self.leaves])
    return(site_pattern)
  def Store_pos_list_for_unique_site_patterns(self):
    self.Set_pre_order_and_post_order()
    self.first_pos_of_unique_site_pattern_to_pos_list = {}
    unique_site_pattern_to_first_pos = {}
    for site in range(self.sequence_length):
      site_pattern = "".join([l.sequence[site] for l in self.leaves])
      # print(f"site pattern for pos{site} is {site_pattern}")
      if site_pattern in unique_site_pattern_to_first_pos.keys():                
        first_pos = unique_site_pattern_to_first_pos[site_pattern] 
        self.first_pos_of_unique_site_pattern_to_pos_list[first_pos].append(site)
      else:
        unique_site_pattern_to_first_pos[site_pattern] = site
        self.first_pos_of_unique_site_pattern_to_pos_list[site] = [site]    
  def Is_informative(self,site_pattern):
    counts = Counter(site_pattern)
    return(sum([count > 1 for count in counts.values()])>1)
  def Get_total_number_of_informative_site_patterns(self):
    total_number_of_informative_sites = 0
    total_number_of_informative_site_patterns = 0
    for site in self.first_pos_of_unique_site_pattern_to_pos_list.keys():     
      site_pattern = "".join([l.sequence[site] for l in self.leaves])  
      if (self.Is_informative(site_pattern)):
        pattern_weight = len(self.first_pos_of_unique_site_pattern_to_pos_list[site])        
        total_number_of_informative_sites += pattern_weight
        total_number_of_informative_site_patterns += 1
    print(f"total number of informative sites is {total_number_of_informative_sites}")
    print(f"total number of unique informative_site_patterns is {total_number_of_informative_site_patterns}")
  def Run_Sankoff(self,site):     
    # for site in range(self.self.sequence_length):
    # Phase 1 compute minimum mutation cost for subtrees
    pattern_weight = len(self.first_pos_of_unique_site_pattern_to_pos_list[site])
    pos_list = self.first_pos_of_unique_site_pattern_to_pos_list[site]
    for parent in self.post_order_list:
      if parent.is_leaf:
        parent.min_cost_subtree = [self.max_penalty] * 4
        dna_p_ind = self.DNA_ind_map[parent.sequence[site]]
        parent.min_cost_subtree[dna_p_ind] = 0        
      else:
        parent.min_cost_subtree = [0] * 4
        for dna_p_ind in range(4):
          for child in parent.children:                        
            min_cost_from_child_subtree = self.max_penalty
            for dna_c_ind in range(4):
              min_cost_from_child_subtree = min(min_cost_from_child_subtree, self.mut_cost_matrix[dna_p_ind][dna_c_ind] + child.min_cost_subtree[dna_c_ind])          
            parent.min_cost_subtree[dna_p_ind] += min_cost_from_child_subtree            
    # Phase 2 compute ancestral sequences
    min_cost_ind = self.root.min_cost_subtree.index(min(self.root.min_cost_subtree))
    self.total_parsimony_score += self.root.min_cost_subtree[min_cost_ind] * pattern_weight    
    dna_root = self.DNA_list[min_cost_ind]
    for pos in self.first_pos_of_unique_site_pattern_to_pos_list[site]:      
      self.root.sequence = self.root.sequence[:pos] + dna_root + self.root.sequence[(pos+1):]
    for child in self.pre_order_list[1:]:
      if not child.is_leaf:
        parent = child.parent
        char_parent = parent.sequence[site]
        dna_p_ind = self.DNA_ind_map[char_parent]
        assert(char_parent != "N")
        min_cost_ind = 0
        min_cost = self.max_penalty
        for dna_c_ind in range(4):
          if self.mut_cost_matrix[dna_p_ind][dna_c_ind] + child.min_cost_subtree[dna_c_ind] < min_cost:
            min_cost = self.mut_cost_matrix[dna_p_ind][dna_c_ind] + child.min_cost_subtree[dna_c_ind]
            min_cost_ind = dna_c_ind
        dna_child = self.DNA_list[min_cost_ind]
        for pos in self.first_pos_of_unique_site_pattern_to_pos_list[site]:
          child.sequence = child.sequence[:pos] + dna_child + child.sequence[(pos+1):]    
  def Run_Sankoff_for_all_sites(self):
    self.Set_pre_order_and_post_order()
    self.Set_null_sequences_in_ancestors()
    for site in self.first_pos_of_unique_site_pattern_to_pos_list.keys():
      site_pattern = "".join([l.sequence[site] for l in self.leaves])
      self.Run_Sankoff(site)
  def Suppress_the_root(self):
    root = self.Get_root()
    self.Set_neighbors(root)
    del self.vertex_map[root.name]
    for v in self.vertex_map.values():
      v.parent = -1
      v.children = []
  def Set_neighbors(self, vertex):
    vertex.neighbors.extend(vertex.children)
    for c in vertex.children:
      if vertex.is_root:
        # c.neighbors.extend(list(set(vertex.children) - set(c)))
        c.neighbors.extend([item for item in vertex.children if item.name != c.name])
        self.Set_neighbors(c)
      else:
        c.neighbors.append(vertex)
        self.Set_neighbors(c)

  def Root_tree_along_branch(self, node1_name, node2_name):
    self.Add_vertex("root")
    self.Set_root()
    self.root.out_degree = 2
    vertices_to_visit = [self.root]
    while len(vertices_to_visit) > 0:
      v = vertices_to_visit.pop()
      if v.is_root:
        v.children = [self.Get_vertex(node1_name), self.Get_vertex(node2_name)]
        for child in v.children:
            child.parent = v

        vertices_to_visit += v.children
      elif v.name == node1_name:
        if v.is_leaf:
          v.neighbors = []
          continue
        else:
          v.children = v.neighbors
          v.children.remove(self.Get_vertex(node2_name))
          for child in v.children:
            child.parent = v
          v.neighbors = []

          vertices_to_visit += v.children
      elif v.name == node2_name:
        if v.is_leaf:
          v.neighbors = []
          continue
        else:
          v.children = v.neighbors
          v.children.remove(self.Get_vertex(node1_name))
          for child in v.children:
            child.parent = v
          v.neighbors = []

          vertices_to_visit += v.children
      else:
        if v.is_leaf:
          v.neighbors = []
          continue
        else:
          v.children = v.neighbors
          v.children.remove(v.parent)
          for child in v.children:
            child.parent = v
          v.neighbors = []

          vertices_to_visit += v.children
    self.total_parsimony_score = 0
    # for v in T.vertex_map.values():
    #   v.neighbors = []
    

  def NNI(self, node1_name, node2_name):
    T1 = copy.deepcopy(self)

    node1_subtrees = T1.Get_vertex(node1_name).neighbors
    # print([v.name for v in node1_subtrees])
    node1_subtrees.remove(T1.Get_vertex(node2_name))
    subtree_1 = node1_subtrees[0]
    subtree_2 = node1_subtrees[1]

    node2_subtrees = T1.Get_vertex(node2_name).neighbors
    node2_subtrees.remove(T1.Get_vertex(node1_name))
    subtree_3 = node2_subtrees[0]
    subtree_4 = node2_subtrees[1]
   
    T1.Get_vertex(node1_name).neighbors = [subtree_1, subtree_3, T1.Get_vertex(node2_name)]
    T1.Get_vertex(node2_name).neighbors = [subtree_4, subtree_2, T1.Get_vertex(node1_name)]
    subtree_2.neighbors.remove(T1.Get_vertex(node1_name))
    subtree_2.neighbors.append(T1.Get_vertex(node2_name))
    subtree_3.neighbors.remove(T1.Get_vertex(node2_name))
    subtree_3.neighbors.append(T1.Get_vertex(node1_name))

    T2 = copy.deepcopy(self)

    node1_subtrees = T2.Get_vertex(node1_name).neighbors
    node1_subtrees.remove(T2.Get_vertex(node2_name))
    subtree_1 = node1_subtrees[0]
    subtree_2 = node1_subtrees[1]

    node2_subtrees = T2.Get_vertex(node2_name).neighbors
    node2_subtrees.remove(T2.Get_vertex(node1_name))
    subtree_3 = node2_subtrees[0]
    subtree_4 = node2_subtrees[1]

    T2.Get_vertex(node1_name).neighbors = [subtree_1, subtree_4, T2.Get_vertex(node2_name)]
    T2.Get_vertex(node2_name).neighbors = [subtree_3, subtree_2, T2.Get_vertex(node1_name)]
    subtree_2.neighbors.remove(T2.Get_vertex(node1_name))
    subtree_2.neighbors.append(T2.Get_vertex(node2_name))
    subtree_4.neighbors.remove(T2.Get_vertex(node2_name))
    subtree_4.neighbors.append(T2.Get_vertex(node1_name))

    return T1, T2


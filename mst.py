# import config
from config import *
#I. Computing pairwise Jukes-Cantor distances between sequences
def jukes_cantor_distance(sequence1, sequence2):
  # print('computing JC distance')
  nucleotide_diff_count = 0 # how many nucleotides they differ
  total_nucleotides_compared = 0.0
  bases = {'A', 'C', 'G', 'T'}
  
  for a, b in zip(sequence1, sequence2):
    if a in bases and b in bases:
      # total_nucleotides_compared += 1
      if a != b: 
        nucleotide_diff_count += 1
      
  nucleotide_diff_percentage = nucleotide_diff_count / len(sequence1)
  # diff_pct_adjusted = min(nucleotide_diff_percentage, 0.74999)
  # jukes_cantor_distance = -0.75 * math.log(1 - min((diff_pct_adjusted*(4.0/3.0)),1.0))
  jukes_cantor_distance = -0.75 * math.log(1 - (nucleotide_diff_percentage*(4.0/3.0)))
  return jukes_cantor_distance

#II. Creating a matrix of pairwise Jukes-Cantor distance
def create_distance_matrix(seq_dict):
  dist_mat_new = np.zeros((len(seq_dict), len(seq_dict)))
  
  for i in range(len(seq_dict)):
    for j in range(i+1, len(seq_dict)):
      dist_mat_new[i][j] = jukes_cantor_distance(seq_dict[list(seq_dict.keys())[i]], seq_dict[list(seq_dict.keys())[j]])
  return dist_mat_new

#III. Creating networkx-based MST with distance matrix
def create_mst_from_dist_mat(seq_dict, distance_matrix):
  seq_id_new = list(seq_dict.keys())
  G_new = nx.Graph()
  G_new.add_nodes_from(seq_id_new)
  for i in range(len(seq_id_new)):
    for j in range(i+1, len(seq_id_new)):
      G_new.add_edge(seq_id_new[i], seq_id_new[j], weight = distance_matrix[i][j])
  MST = nx.minimum_spanning_tree(G_new)
  return MST

from config import *
from tree import *
from mst import *
from infer_local import *
from infer_global import *

def main():

	#Record the starting time
	start_time = time.time()

	#Read the sequences
	sequences = {}
	jc_distances = {}

	fasta_file = open(fasta_file_path,'r')
	seq = ''
	seq_id = ''
	for line in fasta_file:
	  if line.startswith('>'):
	    seq_id = line.strip().split('>')[1]
	    seq = ''
	  else:
	    seq += line.strip()
	  if seq != '' and seq_id != '':
	    sequences[seq_id] = seq

	#Create the pairwise distance matrix
	dist_mat = create_distance_matrix(sequences)

	#Create the networkx-based MST
	MST = create_mst_from_dist_mat(sequences, dist_mat)

	fig, ax = plt.subplots(figsize=(30, 30))
	nx.draw(MST, with_labels=True, font_weight='bold', node_size = 500)
	fig.savefig('MST_starting.png')

	#Create an unrooted global phylogeny tree
	T_global = Tree("global")
	T_global.Read_fasta_file(fasta_file_path)


	test = infer_global_tree(T_global, MST, sequences, W, threshold)
	df.to_csv(f'{len(sequences)}_{threshold}_stats.csv', encoding='utf-8', index=False)
	print("Elapsed time: %s seconds" % (time.time() - start_time))
    
if __name__ == "__main__":

    main()

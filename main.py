# import config
# import tree
# import mst
# import infer_local
# import infer_global
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
	mst_fn = fn + '.MST_starting.png'
	fig.savefig(mst_fn)

	#Create an unrooted global phylogeny tree
	T_global = Tree("global")
	T_global.Read_fasta_file(fasta_file_path)


	T_global = infer_global_tree(T_global, MST, sequences, W, threshold)
	print("Elapsed time: %s seconds" % (time.time() - start_time))

	stat_fn = fn + '.stats.csv'
	df.to_csv(stat_fn, encoding='utf-8', index=False)
	

	seq_df = pd.DataFrame(columns=['Name', "Sequence"])
	i = 1
	for v_name, v in T_global.vertex_map.items():
		seq_df.loc[i] = [v_name, v.sequence]
		i += 1

	seq_fn = fn + '.global_tree_sequences.csv'
	seq_df.to_csv(seq_fn, encoding='utf-8', index=False)

	newick = newick_global_tree(T_global)
	print(newick)
    
if __name__ == "__main__":

    main()

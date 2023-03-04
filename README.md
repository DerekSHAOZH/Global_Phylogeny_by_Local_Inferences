# Global_Phylogeny_by_Local_Inferences
Final Project for FUB 22Win Phylogeny Inference and Application

## .ipynb
*GPLI.ipynb* is the interactive jupyter notebook (supported by Google Colab) which has the same function as the python scripts.

## Prerequisites
* Python (v >= 3.8)

  Some functions in this implementation depend on the insertion-order preservation nature of dict objects which has been declared to be an official part of the Python language spec since Python 3.7.0. 
  
  NetworkX requires Python 3.8, 3.9, or 3.10
  
  You can check your Python version by 
  ```sh
  python3 --version
  ```
* Pandas
  ```sh
  pip install pandas
  ```
  if you are a Conda user,
  ```sh
  conda install -c conda-forge pandas
  ```
  
* Matplotlib
  ```sh
  pip install matplotlib
  ```
  if you are a Conda user,
  ```sh
  conda install -c conda-forge matplotlib
  ```

* Ete3
  ```sh
  pip install ete3
  ```
  if you are a Conda user,
  ```sh
  conda install -c conda-forge ete3
  ```
* Biopython
  ```sh
  pip install biopython
  ```
  if you are a Conda user,
  ```sh
  conda install -c conda-forge biopython
  ```
* Networkx (v=3.0.0)
  ```sh
  pip install networkx
  ```
  
## Installation
```sh
git clone https://github.com/DerekSHAOZH/Global_Phylogeny_by_Local_Inferences/
```

## Execution
* Change dataset path in config.py
  ```python
  #II. Loading the sequence dataset
  fasta_file_path = os.getcwd() + "/Dataset/influenza_98.fasta"  #TODO if change to another dataset
  num_seq = 98    #TODO if change to another dataset
  ```
* Change threshold in config.py
  ```python
  #IV. Setting the threshold for inferring local tree
  threshold = 10  #TODO
  ```  

* Run scripts
  
  Your should run the Python script within this working directory.
  ```sh
  python3 main.py
  ```
## Result
All output files will be directed to the *Result* folder.

For example,

    .
    ├── ...
    ├── Result                                        # files with prefix "98_10" means running with influenza_98 dataset at threshold 10
    │   ├── 98_10.MST_31.png                          # MST graph with 31 nodes
    │   ├── 98_10.MST_65.png                          # MST graph with 65 nodes
    │   ├── 98_10.MST_starting.png                    # Starting MST graph (with 98 nodes)
    │   ├── 98_10.global_tree_sequences.csv           # Sequence of each vertex in global tree
    │   ├── 98_10.log                                 # log file with running time and newick format of global tree
    │   ├── 98_10.stats.csv                           # Statistics per updating iteration
    │   └── 98_10.tmp.tre                             # Temporary file (not important)
    └── ...

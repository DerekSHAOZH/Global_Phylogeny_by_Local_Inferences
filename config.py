#I. Importing required packages
from re import search, compile
import math
import numpy as np
import itertools
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import copy
import os
import ete3
import random
import pandas as pd
import sys
import re
import time
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import *
from io import StringIO

#only in google colab
# from google.colab import drive
# drive.mount('/content/drive')

#II. Loading the sequence dataset
fasta_file_path = "/Dataset/influenza_98.fasta"  #TODO if change to another dataset
num_seq = 98    #TODO if change to another dataset

#III. Defining the mutation cost matrix for maximum parsimony 
W = np.array([[0, 1, 1, 1], 
              [1, 0, 1, 1],
              [1, 1, 0, 1],
              [1, 1, 1, 0]])

#IV. Setting the threshold for inferring local tree
threshold = 10  #TODO

#V. Setting the random seed for reproducibility
random.seed(1234)

#VI. Defining file name
#for example, files with prefix "98_10" means running with influenza_98 dataset at threshold 10
fn = '/Result/' + str(num_seq) + '_' + str(threshold)

#VII. Printing messages to a log file
#not used in .ipynb
orig_stdout = sys.stdout
log_fn = fn + '.log'
sys.stdout = open(log_fn, 'w')
# sys.stdout = orig_stdout

#VIII. Creating a data frame storing stats per iteration
df = pd.DataFrame(columns=['k', 'MP score of local', 'MST nodes', "Global tree vertices", "Global tree edges", "V size"])

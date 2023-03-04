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
from google.colab import drive
drive.mount('/content/drive')

#II. Loading the sequence dataset
fasta_file_path = "/content/drive/MyDrive/WS_22_23/input_files/influenza_98.fasta"  #TODO if change to another dataset
num_seq = 98    #TODO if change to another dataset

#III. Defining the mutation cost matrix for maximum parsimony 
W = np.array([[0, 1, 1, 1], 
              [1, 0, 1, 1],
              [1, 1, 0, 1],
              [1, 1, 1, 0]])

#VII. Setting the threshold for inferring local tree
threshold = 10

#IV. Setting the random seed for reproducibility
random.seed(1234)

#V. Defining file name
#for example, files with prefix "98_10" means running with influenza_98 dataset at threshold 10
fn = '/Result/' + str(num_seq) + '_' + str(threshold)

#VI. Printing messages to a log file
#not used in .ipynb
orig_stdout = sys.stdout
log_fn = fn + '.log'
sys.stdout = open(log_fn, 'w')
# sys.stdout = orig_stdout

#VII. Creating a data frame storing stats per iteration
df = pd.DataFrame(columns=['k', 'MP score of local', 'MST nodes', "Global tree vertices", "Global tree edges", "V size"])

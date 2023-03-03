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

#II. Loading the sequence dataset
fasta_file_path = "/content/drive/MyDrive/WS_22_23/input_files/influenza_98.fasta"

#III. Defining the mutation cost matrix for maximum parsimony 
W = np.array([[0, 1, 1, 1], 
              [1, 0, 1, 1],
              [1, 1, 0, 1],
              [1, 1, 1, 0]])

#IV. Setting the random seed for reproducibility
random.seed(1234)

#V. Printing messages to a log file
orig_stdout = sys.stdout
sys.stdout = open('test.log', 'a')
# sys.stdout = orig_stdout

#VI. Creating a data frame storing stats per iteration
df = pd.DataFrame(columns=['k', 'MP score of local', 'MST nodes', "Global tree vertices", "Global tree edges", "V size"])

#VII. Setting the threshold for inferring local tree
threshold = 10
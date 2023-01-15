import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from protein import Protein
import pandas as pd
import igraph as ig

import networkx as nx

import matplotlib.pyplot as plt

proteins_df = pd.read_csv('data/protein_scores.csv')
interactions_df = pd.read_csv('data/interactions.csv')
interactions_binary = pd.read_csv('data/interactions_binary.csv')
interactions_binary.sort_values(by='interactor_A', inplace=True)

# proteins in the MS data
proteins_1 = proteins_df.iloc[:, 0]
# proteins in the interactions data
proteins_2 = interactions_df.iloc[:, 0]
all_proteins = pd.concat([proteins_1, proteins_2], axis=0)
all_proteins.reset_index(drop=True)


proteins = {}


# def append_to_dict(x):
#     proteins[x['protein_accession']] = Protein(x['protein_accession'], x['protein_score'])

# proteins_df.apply(append_to_dict, axis=1)

def append_to_dict(x):
    proteins[x['interactor_A']] = Protein(accession=x['interactor_A'], interactions=x['interact_with'], score=0)

# create protein objects from pandas dataframe
interactions_df.apply(append_to_dict, axis=1)

# update the protein score information from the MS data
def update_score(x):
    # check if the protein has interaction information
    if x['protein_accession'] in proteins.keys():
        proteins[x['protein_accession']].score = x['protein_score']
    # if there is no interaction information
    else:
        proteins[x['protein_accession']] = Protein(accession=x['protein_accession'], interactions=None, score=x['protein_score'])

proteins_df.apply(update_score, axis=1)

print(proteins['P07355'])


edges = []
def create_edge_list_from_row(x):
    edges.append((x['interactor_A'], x['interactor_B']))

interactions_binary.apply(create_edge_list_from_row, axis=1)
print(edges)


G = nx.Graph(edges)
nx.draw(G, node_size=10)
plt.show()
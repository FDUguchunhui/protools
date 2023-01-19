import torch
from torch_geometric.data import InMemoryDataset, Data
import pandas as pd

class ProteinDataSet(InMemoryDataset):
    def __init__(self, root, transform=None, pre_transform=None, pre_filter=None):
        super().__init__(root, transform, pre_transform, pre_filter)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return []

    @property
    def processed_file_names(self):
        return ['data.pt']

    def process(self):

        proteins = pd.read_csv(self.root + 'protein-dataset.csv')





        #
        # if self.pre_filter is not None:
        #     data_list = [data for data in data_list if self.pre_filter(data)]
        #
        # if self.pre_transform is not None:
        #     data_list = [self.pre_transform(data) for data in data_list]





        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])
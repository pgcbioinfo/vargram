"""Main module for data wrangling."""

from ._nextclade import nextclade
from . import _nextclade_utils
import pandas as pd
import os

def read_table(table_object, nextclade_file=False):
    if isinstance(table_object, pd.DataFrame):
        table = table_object
    else:
        ext = os.path.splitext(table_object)[1]
        if nextclade_file:
            csv_delim = ';'
        else:
            csv_delim = ','

        match ext:
            case '.csv':
                delimiter = csv_delim
            case '.tsv':
                delimiter = '\t'
            case _:
                raise ValueError(f"Unrecognized file path. Expecting .csv or .tsv file but got {ext}")
        
        table = pd.read_csv(table_object, delimiter=delimiter) 
    return table

class Wrangler():

    def __init__(self, wrangler_kwargs):
        """Determine which plot to process data for."""
        self.wrangler_kwargs = wrangler_kwargs
        self.plot = wrangler_kwargs['plot']
        self.format = wrangler_kwargs.get('format')
        del wrangler_kwargs['plot']
        if self.format is not None:
            del wrangler_kwargs['format']
        self.user_input = wrangler_kwargs

        match self.plot:
            case 'Profile':
                self._profile()
        
        if self.data.empty:
            raise ValueError("Wrangled data is empty.")
    
    def get_wrangled_data(self):
        return self.data, self.format
            
    def _profile(self):
        """Choosing appropriate data wrangling method for Profile()."""
        profile_formats = ['nextclade_fasta', 'nextclade_delimited', 'vargram', 'delimited', '_test']

        # Assigning default format values
        if self.format is None:
            if 'seq' in self.user_input.keys():
                self.format = 'nextclade_fasta'
            elif 'data' in self.user_input.keys():
                self.format = 'nextclade_delimited'
        elif self.format not in profile_formats:
            raise ValueError(f"Unrecognized format: {self.format}.")

        # Wrangling data
        match self.format:
            case 'nextclade_fasta':
                nextclade_kwargs = {key: self.user_input[key] for key in ['seq', 'ref', 'gene'] if key in self.user_input.keys()}
                read_data = nextclade(**nextclade_kwargs)
                self.data = _nextclade_utils.process_nextclade(read_data)
            case 'nextclade_delimited':
                tabular_data = self.user_input['data']
                read_data = read_table(tabular_data, nextclade_file=True)
                if 'batch' not in read_data.columns:
                    read_data.insert(0, 'batch', 'my_batch')
                read_data.sort_values(by=['batch', 'seqName'], inplace=True)
                read_data.reset_index(drop=True, inplace=True)
                self.data = _nextclade_utils.process_nextclade(read_data)
            case _:
                tabular_data = self.user_input['data']
                read_data = read_table(tabular_data)
                self.data = read_table(tabular_data, nextclade_file=False)

        # Adding metadata if available
        if 'meta' in self.user_input.keys():
            # Getting columns to merge on 
            if 'join' not in self.user_input.keys():
                missing_join_error=("Missing 'join' argument. Provide column name(s) " 
                                    "to (outer) merge data and metadata on.")
                raise ValueError(missing_join_error)
            elif isinstance(self.user_input['join'], str): 
                join = [self.user_input['join']]
            else:
                join = self.user_input['join']
            metadata = self.user_input['meta']
            if len(join) != 1:
                metadata.rename(columns={join[1]:join[0]}, inplace=True)
                self.data = pd.merge(self.data, metadata, how='outer', on=join[0])
            elif self.format == 'nextclade_fasta' or self.format == 'nextclade_delimited': 
                # If 'seq' is provided, automatically join on nextclade_sequence_name
                nextclade_seqname = 'seqName'
                metadata.rename(columns={join[0]:nextclade_seqname}, inplace=True)
                self.data = pd.merge(self.data, metadata, 
                                      how='outer', on=nextclade_seqname)
            else: 
                self.data = pd.merge(self.data, metadata, how='outer', on=join)
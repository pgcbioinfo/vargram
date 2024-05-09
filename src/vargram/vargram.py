# """Main module for generating VARGRAM figures and summary statistics."""

from .methods import stats_module, stats_wrap, bar_module, bar_wrap, methods_utils
from .nextread.nextread import nextread
from .plots.bar import Bar

import matplotlib.pyplot as plt
import pandas as pd

class vargram:

    def __init__(self, **vargram_kwargs):
        """Initializes the VARGRAM object, getting data and joining it with metadata (if provided).

        Args:
            seq (string): FASTA file path of the sequences.
            ref (string): FASTA file path of the reference sequence.
            gene (string): GFF3 file path of the genome annotation.
            data (string or pandas.DataFrame): The data to be plotted.
            metadata (string or pandas.DataFrame): The metadata to be joined with the data.
            join (string or list): Provides the column to do an outer merge on. 
                If data and metadata do not share a column, the two column names which share the values of interest may be provided.

        Raises:
            ValueError: The resulting combined data and metadata is empty.
        """

        # Defining initial values
        self.__initialize_variables__()

        # Getting data
        if 'seq' in vargram_kwargs.keys(): # Get Nextclade output
            nextread_kwargs = {key: vargram_kwargs[key] for key in ['seq', 'ref', 'gene'] if key in vargram_kwargs.keys()}
            self._data = nextread(**nextread_kwargs)

            # Defining relevant Nextclade analysis column names
            self._nextclade_seqname = 'seqName'
            self._aasub = 'aaSubstitutions'
            self._aadel = 'aaDeletions'
            self._aains = 'aaInsertions'
            self._nextclade_seqname = 'seqName'
        elif 'data' not in vargram_kwargs.keys():
            raise ValueError("Missing data. Either provide the FASTA files with the 'seq' and 'ref' arguments or provide a dataframe/path with 'data'.")
        else:
            if isinstance(vargram_kwargs['data'], str): # data is a path
                self._data = pd.read_csv(vargram_kwargs['data'])
            else:
                self._datdata = vargram_kwargs['data'] # data is a dataframe
        
        # Getting metadata and joining it with data
        if 'meta' in vargram_kwargs.keys():

            # Getting columns to merge on 
            if 'join' not in vargram_kwargs.keys():
                raise ValueError("Missing 'join' argument. Provide column name(s) to (outer) merge data and metadata on.")
            elif isinstance(vargram_kwargs['join'], str): 
                join = [vargram_kwargs['join']]
            else:
                join = vargram_kwargs['join']
            
            metadata = vargram_kwargs['meta']
            if len(join) != 1:
                metadata.rename(columns={join[1]: join[0]}, inplace=True)
                self._data = pd.merge(self._data, metadata, how='outer', on=join[0])
            elif 'seq' in vargram_kwargs.keys(): # If 'seq' is provided, automatically join on nextclade_sequence_name
                metadata.rename(columns={join[0]: self._nextclade_seqname}, inplace=True)
                self._data = pd.merge(self._data, metadata, how='outer', on=self._nextclade_seqname)
            else: 
                self._data = pd.merge(self._data, metadata, how='outer', on=join)                

        if self._data.empty:
            raise ValueError("Data is empty.")
    
    def __initialize_variables__(self):
        self._methods_called = []
        self._methods_kwargs = []
        self._plots = ['bar'] # These are methods that build the plot
        self._terminals = ['show', 'savefig']
        self._generate_plot = False

        # show()
        self._shown = False

        # savefig()
        self._saved = False
    
    def _generate(self):
        
        # Getting group of methods called since last called plot
        latest_method_calls = self._methods_called[self._recent_plot_index:]

        # Determining whether to generate plots or not
        if '_show' in latest_method_calls and '_savefig' in latest_method_calls:
            getattr(self, self._methods_called[-1])(**self._methods_kwargs[-1])
            return

        # Generating figure based on most recent plot called
        for i in range(self._recent_plot_index, len(self._methods_called)):
            # Running each method since the most recent plot called
            getattr(self, self._methods_called[i])(**self._methods_kwargs[i])

        
    
    def _show(self, empty_string=''): # The unused empty string argument is so as to be able to maintain length of methods and methods_kwargs the same

        print("Showed figure.")#plt.show()
    
    def _savefig(self, **_savefig_kwargs):

        print("Saved figure.")#plt.savefig(**_savefig_kwargs, bbox_inches='tight_layout')
    
    def _bar(self, **_bar_kwargs):

        Bar(self._data, **_bar_kwargs)

    def show(self): 
        self._shown = True
        self._methods_called.append('_show')
        self._methods_kwargs.append({'empty_string':''})

        self._generate()

    def savefig(self, **savefig_kwargs):
        self._saved = True
        self._methods_called.append('_savefig')
        self._methods_kwargs.append(savefig_kwargs)

        self._generate()

    def bar(self, **bar_kwargs):
        self._methods_called.append('_bar')
        self._methods_kwargs.append(bar_kwargs)
        self._recent_plot_index = len(self._methods_called) - 1
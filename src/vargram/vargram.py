# """Main module for generating VARGRAM figures and summary statistics."""

from .utils import methods_utils
from .nextread.nextread import nextread
from .nextread import nextread_utils
from .plots.bar import Bar

import pandas as pd
import os

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
        self._initialize_variables()

        # Getting data
        if 'seq' in vargram_kwargs.keys(): # Get Nextclade output
            # Defining relevant Nextclade analysis column names
            self._nextclade_seqname = 'seqName'
            self._aasub = 'aaSubstitutions'
            self._aadel = 'aaDeletions'
            self._aains = 'aaInsertions'
            self._nextclade_seqname = 'seqName'

            nextread_kwargs = {key: vargram_kwargs[key] for key in ['seq', 'ref', 'gene'] if key in vargram_kwargs.keys()}
            self._data = nextread(**nextread_kwargs)
            self._data = nextread_utils.process_nextread(self._data)
            self._nextread_called = True
            
        elif 'data' not in vargram_kwargs.keys():
            raise ValueError("Missing data. Either provide the FASTA files with the 'seq' and 'ref' arguments or provide a dataframe/path with 'data'.")
        else:
            if isinstance(vargram_kwargs['data'], str): # data is a path

                ext = os.path.splitext(vargram_kwargs['data'])[1]
                if vargram_kwargs['nextclade'] == True:
                    delimiter = ';'
                elif ext == '.csv':
                    delimiter = ','
                elif ext == '.tsv':
                    delimiter = '\t'

                read_data = pd.read_csv(vargram_kwargs['data'], delimiter=delimiter)
            else:
                read_data = vargram_kwargs['data'] # data is a dataframe
            
            # Processing if provided CSV is Nextclade output
            if vargram_kwargs['nextclade'] == True:
                if 'batch' not in read_data.columns:
                    read_data.insert(0, 'batch', 'my_batch')
                self._data = nextread_utils.process_nextread(read_data)
            else:
                self._data = read_data
        
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
    
    def _initialize_variables(self):
        """Sets initial values of attributes.
        """
        self._methods_called = []
        self._methods_kwargs = []
        self._plots = ['bar'] # These are methods that build the plot
        self._terminals = ['show', 'save', 'stat']
        self._generate_plot = False
        self._plots_indices = []

        # key()
        self._clean_keys()

        # show()
        self._shown = False

        # save()
        self._saved = False
    
    def _clean_keys(self):

        self._master_key_data = pd.DataFrame()
        self._nkeys = 0
        self._key_labels = []
        self._key_colors = []


    
    def _generate(self):
        """Runs all called methods in correct order"""
        
        # Getting group of methods called since last called plot
        latest_method_calls = self._methods_called[self._latest_plot_index:]
        latest_method_kwargs = self._methods_kwargs[self._latest_plot_index:]

        # Determining whether to generate plots or not
        if not self._generate_plot:
            if latest_method_calls[-1] == '_stat':
                return
            else:
                getattr(self, latest_method_calls[-1])(**self._methods_kwargs[-1])
                return
        
        # Creating plot object instance
        plot_class = latest_method_calls[0][1:].title() # This removes the initial underscore and capitalizes the first letter
        plot_object = globals()[plot_class]
        self._plot_instance = plot_object(self._data)
        
        # Rearranging so that auxiliary methods are run before plot and save/show methods
        latest_method_calls.append(latest_method_calls[-1])
        latest_method_calls[-2] = latest_method_calls[0]
        latest_method_calls.pop(0)
        latest_method_kwargs.append(latest_method_kwargs[-1])
        latest_method_kwargs[-2] = latest_method_kwargs[0]
        latest_method_kwargs.pop(0)

        # Generating figure based on most recent plot called
        run_key = True
        for i, method in enumerate(latest_method_calls):
            # Running each method since the most recent plot called
            if run_key or method != '_key': # Ensuring that the internal key method is only called once
                if method == '_key':
                    run_key = False
                
                if method == '_bar':
                    self._plot_data = getattr(self, method)(**latest_method_kwargs[i])
                else:
                    getattr(self, method)(**latest_method_kwargs[i])
            else:
                continue
        
        self._generate_plot = False
    
    def _show(self, empty_string=''): 
        """Show generated figure"""

        getattr(self._plot_instance, 'show')()
    
    def _save(self, **_save_kwargs):
        """Save generated figure or data"""

        getattr(self._plot_instance, 'save')(**_save_kwargs)
    
    def _stat(self, **_stat_kwargs):
        
        pass
    
    def _bar(self, **_bar_kwargs):
        """Process data for plotting"""
        
        return getattr(self._plot_instance, 'process')(**_bar_kwargs)

    def _key(self, **_key_kwargs):
        """Process key data for plotting"""

        getattr(self._plot_instance, 'key')(key_data = self._master_key_data, 
                                            key_labels = self._key_labels, 
                                            key_colors = self._key_colors)
    
    def _struct(self, **_struct_kwargs):
        """Modify aesthetic attributes"""
        
        getattr(self._plot_instance, 'struct_method')(**_struct_kwargs)

    def _aes(self, **_aes_kwargs):
        """Modify aesthetic attributes"""
        
        getattr(self._plot_instance, 'aes')(**_aes_kwargs)
    
    def stat(self, **stat_kwargs):
        """Wrapper for stat method"""

        self._methods_called.append('_stat')
        self._methods_kwargs.append({'empty_string':''}) # The unused empty string argument is so as to be able to maintain length of methods and methods_kwargs the same

        self._generate()

        self._plot_data.sort_values(by=['gene', 'mutation'], inplace=True)
        self._plot_data.reset_index(drop=True, inplace=True)

        return self._plot_data
        
    def show(self): 
        """Wrapper for show method"""

        self._shown = True
        self._methods_called.append('_show')
        self._methods_kwargs.append({'empty_string':''}) # The unused empty string argument is so as to be able to maintain length of methods and methods_kwargs the same

        self._generate()

    def save(self, fname, **save_kwargs):
        """Wrapper for save method"""

        self._saved = True
        save_kwargs['fname'] = fname
        self._methods_called.append('_save')
        self._methods_kwargs.append(save_kwargs)

        self._generate()

    def bar(self, **bar_kwargs):
        """Captures bar method arguments"""
        
        self._methods_called.append('_bar')
        self._methods_kwargs.append(bar_kwargs)
        self._latest_plot_index = len(self._methods_called) - 1
        self._generate_plot = True
        self._clean_keys()

    def key(self, key_data, **key_kwargs):
        """Joins all key data"""

        self._methods_called.append('_key')
        self._methods_kwargs.append({'empty_string':''}) # The unused empty string argument is so as to be able to maintain length of methods and methods_kwargs the same
        self._nkeys += 1

        # Reading data
        if isinstance(key_data, str):
            key_read = methods_utils.read_comma_or_tab(key_data)

        # Getting x to read
        if 'x' in key_kwargs.keys():
            key_x = key_kwargs['x']
        else:
            key_x = 'mutation'
        
        # Getting group to read
        if 'group' in key_kwargs.keys():
            key_group = key_kwargs['group']
        else:
            key_group = 'gene'
        
        # Just keeping the 'x' and 'group' columns
        key_df = key_read[[key_group, key_x]].copy()

        # Getting name of lineage
        if 'label' in key_kwargs:
            key_label = key_kwargs['label']
        elif isinstance(key_data, str):
            key_label = os.path.basename(key_data)
            key_label = os.path.splitext(key_label)[0]
        else:
            key_label = f'key_{self._nkeys}'
        self._key_labels.append(key_label)

        # Getting color of lineage
        if 'color' in key_kwargs:
            color = key_kwargs['color']
        else:
            color = '#5E5E5E'
        self._key_colors.append(color)

        # Gathering in one master dataframe
        if self._nkeys > 1:
            key_df[key_label] = 1
            self._master_key_data = pd.merge(self._master_key_data, key_df, on=[key_group, key_x], how='outer')
            self._master_key_data.fillna(0, inplace=True)
            self._master_key_data.reset_index(drop=True, inplace=True)
        else:
            self._master_key_data = pd.DataFrame()
            self._master_key_data[key_group] = key_df[key_group]
            self._master_key_data[key_x] = key_df[key_x]
            self._master_key_data[key_label] = 1

    def struct(self, struct_arg):
        self._methods_called.append('_struct')
        self._methods_kwargs.append({'struct_key':struct_arg})

    def aes(self, **aes_kwargs):
        """Captures modification of aesthetic attributes"""
        
        self._methods_called.append('_aes')
        self._methods_kwargs.append(aes_kwargs)        
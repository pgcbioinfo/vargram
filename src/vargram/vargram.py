"""Main module for generating VARGRAM figures and data."""

from .wranglers._wrangler import Wrangler
from .plots._profile import Profile
import pandas as pd
import os


class vargram:

    def __init__(self, **vargram_kwargs):
        """Initializes the VARGRAM object, getting data 
        and joining it with metadata (if provided).

        Parameters
        ----------
        seq : str
            FASTA file path of the sequences.
        ref : str
            FASTA file path of the reference sequence.
        gene : str
            GFF3 file path of the genome annotation.
        data : str or pandas.DataFrame
            The data to be plotted.
        metadata : string or pandas.DataFrame
            The metadata to be joined with the data.
        join : str or list
            Provides the column to do an outer merge on. 
            If data and metadata do not share a column, 
            the two column names which share the values of interest may be provided.
        
        Returns
        -------
        None

        Raises
        ------
        ValueError
            If neither 'seq' nor 'data' is provided.
            If 'metadata' is provided but 'join' is not.
            If the resulting '_data' attribute is empty.

        """
        # Defining initial values
        self._initialize_variables()
        self._vargram_kwargs = vargram_kwargs
    
    def _initialize_variables(self):
        """Sets initial values of attributes."""
        self._methods_called = [] # List of all methods called
        self._methods_kwargs = [] # List of the arguments of the corresponding methods called
        self._plots = ['profile'] # These are the 'plot methods'
        self._terminals = ['show', 'save', 'stat'] # These are the 'terminal methods'
        self._generate_plot = False
        self._plots_indices = []
        self._clean_keys()
        self._shown = False
        self._saved = False
    
    def _clean_keys(self):
        """Flushes key variables clean."""
        self._master_key_data = pd.DataFrame()
        self._nkeys = 0
        self._key_labels = []
        self._key_colors = []
    
    def _generate(self):
        """Runs all called methods in correct order"""

        # Ensure aes method is called
        if '_aes' not in self._methods_called[self._latest_plot_index:]:
            self.aes() # Call aes to get default arguments
            # Move the aes method name and arguments next to the plot method called
            self._methods_called.insert(self._latest_plot_index + 1, '_aes')
            self._methods_kwargs.insert(self._latest_plot_index + 1, self._methods_kwargs[-1])
            self._methods_called.pop()
            self._methods_kwargs.pop()
        
        # Getting group of methods and corresponding arguments called since last called plot
        latest_method_calls = self._methods_called[self._latest_plot_index:]
        latest_method_kwargs = self._methods_kwargs[self._latest_plot_index:]

        # Determining whether to generate plots or not
        if not self._generate_plot:
            getattr(self, latest_method_calls[-1])(**self._methods_kwargs[-1])
            return
        
        # Creating plot object instance
        plot_class = latest_method_calls[0][1:].title() 
        plot_object = globals()[plot_class]
        self._vargram_kwargs['plot'] = plot_class
        wrangled_data = Wrangler(self._vargram_kwargs).get_wrangled_data()
        self._plot_instance = plot_object(wrangled_data)

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
                if method == '_profile':
                    self._plot_data = getattr(self, method)(**latest_method_kwargs[i])
                else:
                    getattr(self, method)(**latest_method_kwargs[i])
            else:
                continue
        self._generate_plot = False
    
    def _show(self, **_show_kwargs): 
        """Show generated figure"""
        getattr(self._plot_instance, 'show')()
    
    def _save(self, **_save_kwargs):
        """Save generated figure or data"""
        getattr(self._plot_instance, 'save')(**_save_kwargs)
    
    def _stat(self, **_stat_kwargs):
        """Kept for consistency with _generate() flow."""
        pass
    
    def _profile(self, **_profile_kwargs):
        """Process data for plotting"""
        return getattr(self._plot_instance, 'process')(**_profile_kwargs)

    def _key(self, **_key_kwargs):
        """Process key data for plotting"""
        getattr(self._plot_instance, 'key')(key_data=self._master_key_data, 
                                            key_labels=self._key_labels, 
                                            key_colors=self._key_colors)
    
    def _struct(self, **_struct_kwargs):
        """Modify aesthetic attributes"""
        getattr(self._plot_instance, 'struct_method')(**_struct_kwargs)

    def _aes(self, **_aes_kwargs):
        """Modify aesthetic attributes"""
        getattr(self._plot_instance, 'aes')(**_aes_kwargs)
    
    def stat(self):
        """Wrapper for stat method
        
        Returns
        -------
        None

        """
        self._methods_called.append('_stat')
        self._methods_kwargs.append({'empty_string':''}) 
        # The unused empty string argument is so as to be able to maintain 
        # length of methods and methods_kwargs the same
        self._generate()
        return self._plot_data
        
    def show(self): 
        """Wrapper for show method"""
        self._shown = True
        self._methods_called.append('_show')
        self._methods_kwargs.append({'empty_string':''}) 
        # The unused empty string argument is so as to be able to maintain
        # length of methods and methods_kwargs the same
        self._generate()

    def save(self, fname, **save_kwargs):
        """Wrapper for save method"""
        self._saved = True
        save_kwargs['fname'] = fname
        self._methods_called.append('_save')
        self._methods_kwargs.append(save_kwargs)
        self._generate()

    def profile(self,
                threshold=10, 
                x='mutation', 
                y='',
                ytype='', 
                group='gene', 
                stack='batch',
                ):
        """Captures profile method arguments
        
        Parameters
        ----------
        x : str, default:'mutation'
            The x-axis variable (column name of data).
        y : str
            The y-axis variable (column name of data).
        ytype : {'counts', 'weights'}
            The type of calculation for the y-axis.
            'counts' indicate raw counts.
            'weights' indicate normalization by total count.
        group : str, default:'gene'
            The column name for the groups.
        stack : str, default:'batch'
            The column name for the stacks.
        threshold : int, default:50
            The minimum number of occurences of an x value to be included.

        Returns
        -------
        None

        """
        profile_kwargs = locals()
        del profile_kwargs['self']
        self._methods_called.append('_profile')
        self._methods_kwargs.append(profile_kwargs)
        self._latest_plot_index = len(self._methods_called) - 1
        self._generate_plot = True
        self._clean_keys()

    def key(self, key_data, x='mutation', group='gene', label='', color = '#5E5E5E'):
        """Joins all key data.
        
        Parameters
        ----------
        key_data : str
            File path of key lineage
        x : str, default:'mutation'
            The x-axis variable (column name of data).
        group : str, default:'gene'
            The column name for the groups.
        label : str
            The name of the key lineage.
        color : str
            The color of key mutations on the generated heatmap.
        
        Returns
        -------
        None

        """
        self._methods_called.append('_key')
        self._methods_kwargs.append({'empty_string':''}) 
        # The unused empty string argument is so as to be able to maintain 
        # length of methods and methods_kwargs the same

        self._nkeys += 1
        # Reading data
        if isinstance(key_data, str):
            key_read = self._read_comma_or_tab(key_data)
        if isinstance(key_data, pd.DataFrame):
            key_read = key_data.copy()

        # Just keeping the 'x' and 'group' columns
        key_df = key_read[[group, x]].copy()

        # Getting name of key
        if label == '' and isinstance(key_data, str):
            label = os.path.basename(key_data)
            label = os.path.splitext(label)[0]
        self._key_labels.append(label)

        # Getting color of key
        self._key_colors.append(color)

        # Gathering in one master dataframe
        if self._nkeys > 1:
            key_df[label] = 1
            self._master_key_data = pd.merge(self._master_key_data, key_df, 
                                             on=[group, x], how='outer')
            self._master_key_data.fillna(0, inplace=True)
            self._master_key_data.reset_index(drop=True, inplace=True)
        else:
            self._master_key_data = pd.DataFrame()
            self._master_key_data[group] = key_df[group]
            self._master_key_data[x] = key_df[x]
            self._master_key_data[label] = 1

    def struct(self, struct_arg):
        """Retains only groups that are provided.
        
        Parameters
        ----------
        struct_arg : str
            A formatted string indicating the groups to retain per row.
            Commas separate group names and forward slashes separate rows.

        Returns
        -------
        None

        """
        self._methods_called.append('_struct')
        self._methods_kwargs.append({'struct_key':struct_arg})

    def aes(self,
            order=False, 
            flat=False,
            stack_label=[], 
            stack_title='', 
            stack_color='', 
            group_title='', 
            legtitle_fontsize='large', 
            legentry_fontsize='large',
            xticks_fontsize=6, 
            xticks_rotation=90,
            yticks_fontsize=10,
            ylabel='', 
            ylabel_fontsize='large', 
            group_fontsize='large',
            key_fontsize='medium',
            aspect=None
            ):
        """Captures modification of aesthetic attributes.
        
        Parameters
        ----------
        order : bool, default:False
            Determines whether genes will be ordered according to start position.
        flat : bool, default:False
            Determines whether genes will be shown on a horizontal layout.
        stack_label : str
            Name(s) of the stack(s).
        stack_title : str
            Legend title for the stack(s).
        stack_color : str or list
            Color(s) of the stack(s).
        group_title : str
            Legend title for the groups.
        legtitle_fontsize : str or float
            The fontsize of the legend titles.
        legentry_fontsize : str or float
            The fontsize of the legend entries.
        xticks_fontsize : float
            The fontsize of the x-axis ticks.
        xticks_rotation : float
            The degree of rotation of the x-axis ticks.
        yticks_fontsize : float
            The fontsize of the y-axis ticks.
        ylabel : str
            The y-axis label.
        ylabel_fontsize : str or float
            The fontsize of the y-axis label.
        group_fontsize : str or float
            The fontsize of the group label.
        key_fontsize : str or float
            The fontsize of the heatmap (key) row labels.
        aspect : float
            The aspect ratio of the figure.

        Returns
        -------
        None

        """
        aes_kwargs = locals()
        del aes_kwargs['self']
        self._methods_called.append('_aes')
        self._methods_kwargs.append(aes_kwargs)

    def _read_comma_or_tab(self, file):
        """Reads a CSV or TSV file.

        Parameters
        ----------
        file : str
            The file path to be read.
        
        Returns
        -------
        pandas.DataFrame
            The DataFrame of the file.
        
        Raises
        ------
        ValueError
            If the given file is not a CSV or TSV file.

        """
        root, extension = os.path.splitext(file)
        if extension == '.csv':
            return pd.read_csv(file)
        elif extension == '.tsv':
            return pd.read_csv(file, delimiter='\t')
        else:
            raise ValueError("Incorrect file format. Expecting CSV or TSV file.")        
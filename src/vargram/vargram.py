# """Main module for generating VARGRAM figures and summary statistics."""

from .methods import stats_module, stats_wrap, bar_module, bar_wrap, methods_utils

import matplotlib.pyplot as plt
import matplotlib.gridspec as mg
import matplotlib.colors as mc
import pandas as pd
import numpy as np
from collections import Counter
import os

class vargram:

    def __init__(self, nextread_out, threshold = 50, yaxis = []):

        if nextread_out.empty:
            raise ValueError("Analysis DataFrame is empty.")

        self._nextread = nextread_out #[!]Determine if this has use in the future[!]
        self._batch_counts = stats_module.count_per_batch(self._nextread, threshold = threshold)
        self._batches = self._batch_counts['batch'].unique()
        if len(yaxis) == 0:
            if len(self._batches) == 1:
                self._yaxis = 'Instances'
            else:
                self._yaxis = 'Weights'
        else:
            self._yaxis = yaxis
        self._summary_counts = stats_wrap.pivot_batch_counts(self._batch_counts, self._batches, self._yaxis)

        if self._yaxis != 'Instances': # Stats for saving should only be in terms of raw counts
            self._summary_for_saving = stats_wrap.pivot_batch_counts(self._batch_counts, self._batches, yaxis='Instances')

        self._gene_title_axes = []
        self._gene_mutation_axes = []
        self._gene_key_axes = []

        self._show_called = False
        self._bar_called = False
        self._key_called = False
        self._savefig_called = False
        self._stats_called = False
        self._key_labels = []
        self._key_paths = []
        self._key_colors = []
        self._struct = ''
    
    def show(self, suppress = False):
        """Main method of VARGRAM. If this is not called, other methods will not run."""

        self._show_called = True

        # Checking usage
        #if not methods_utils.usage_checker(self._bar_called): #[!]Expand this in the future[!]
        #   raise methods_utils.UsageError("Generate barplot using bar().")        

        # Creating matplotlib Figure object
        self._fig = plt.figure(figsize=[8, 4.8])

        # Creating bar plot
        if self._bar_called:
            self.bar(**self._bar_kwargs)

        # Saving the plot
        if self._savefig_called:
            self.savefig(*self._savefig_args, **self._savefig_kwargs)

        # Saving statistics
        if self._stats_called:
            self.stats(*self._stats_args)

        if not suppress:
            plt.tight_layout()
            plt.show()

    def bar(self, **kwargs):
        """Method that generates the VARGRAM bar plot.

        Usage:
            vg.bar(color = ['white', 'red', 'orange'],
               label = ['Batch 1', 'Batch 2'])
        
        Args:
            color (string or list), optional: List of colors to be used for different batches.
            label (string or list), optional: List of batch names.
        
        Returns:
            None.
        """

        self._bar_called = True
        self._bar_kwargs = kwargs
        
        if self._show_called == False: 
            return None
        
        # Getting structure of the mutation profile grid
        gene_counts = stats_module.count_per_gene(self._keyful_counts)
        if len(self._struct) == 0:
            self._struct = bar_module.build_struct(gene_counts)            
            
        #  Making the barplot
            # Getting bar arguments
        bar_keys = self._bar_kwargs.keys()

            # Getting user-provided colors
        if 'color' in bar_keys:

            # Check if the length of colors matches the number of batches
            if len(self._bar_kwargs['color']) == len(self._batches):
                self._bar_colors = self._bar_kwargs['color']
            else:
                raise ValueError("No. of colors ({}) provided does not match no. of batches ({}).".format(len(self._bar_kwargs['color']), len(self._batches)))
            # Creating default bar colors
        else:
            # Define the colormap
            cmap_name = 'viridis'  
            cmap = plt.get_cmap(cmap_name)

            # Number of colors you want in your list
            num_colors = len(self._batches)

            # Create the ListedColormap
            listed_cmap = cmap(np.linspace(0, 1, num_colors))

            # Convert the colors to hex format (optional)
            self._bar_colors = [mc.to_hex(color) for color in listed_cmap]
        
            # Getting user-provided labels
        if 'label' in bar_keys:
            
            if len(self._bar_kwargs['label']) == len(self._batches):
                self._bar_labels = self._bar_kwargs['label']
            else:
                raise ValueError("No. of labels ({}) provided does not match no. of batches ({}).".format(len(self._bar_kwargs['label']), ))
            
            # Getting user-provided legend title
        if 'legend_title' in bar_keys:
            legend_title =  self._bar_kwargs['legend_title']
        else:
            legend_title = 'Batch'

        # Creating bar grids
        grid_width_counts = stats_module.count_per_gene(self._keyful_counts)
        label_grid, legend_grid, self._gene_title_axes, self._gene_mutation_axes, self._gene_key_axes  = bar_wrap.build_bar_grid(self._struct, grid_width_counts, self._key_called)
        
        # Building mutation profile
        gene_labels = []
        bar_wrap.build_mutation_profile(self._gene_title_axes, self._gene_mutation_axes, self._gene_key_axes,
                                        self._keyful_counts, 
                                        self._struct, 
                                        self._fig, 
                                        self._key_called, self._key_labels, self._key_colors, 
                                        self._batches, self._bar_colors,
                                        gene_labels)
        
        # Creating figure y-axis label
        bar_wrap.build_yaxis_label(self._yaxis, label_grid)

        # Creating figure legend
        bar_wrap.build_legend(legend_grid, self._batches, self._bar_colors, legend_title, gene_labels)
                    
    def key(self, path, **kwargs):

        self._key_called = True
        self._key_paths.append(path)

        if 'label' in kwargs.keys():
            key_label = kwargs['label']
            
        else:
            key_label = os.path.basename(path)
            key_label = os.path.splitext(key_label)[0]

        self._key_labels.append(key_label)

        if 'color' in kwargs.keys():
            self._key_colors.append(kwargs['color'])
        else:
            self._key_colors.append('#5E5E5E')

        if len(self._key_paths) == 1: # Only make a copy once a key lineage is called
            self._keyful_counts = self._summary_counts.copy()
            if self._yaxis != 'Instances': # Making a copy of raw batch counts
                self._keyful_for_saving = self._summary_for_saving.copy()

        current_key = pd.read_csv(path)
        self._keyful_counts = stats_wrap.add_key_mutations(self._keyful_counts, current_key, key_label)
        if self._yaxis != 'Instances': # Making a copy of raw batch counts
            self._keyful_for_saving = stats_wrap.add_key_mutations(self._keyful_for_saving, current_key, key_label)

        if self._show_called == False: 
            return None
        
    def struct(self, struct):

        self._struct = struct
    
    def stats(self, *args):

        self._stats_called = True
        self._stats_args = args

        if self._show_called == False: 
            return None

        if self._yaxis == 'Instances':
            self._keyful_counts.sort_values(by='gene').to_csv(*self._stats_args, index = False)
        else:
            self._keyful_for_saving.sort_values(by='gene').to_csv(*self._stats_args, index = False)

        return None

    def savefig(self, *args, **kwargs):

        self._savefig_called = True
        self._savefig_kwargs = kwargs
        self._savefig_args = args

        if self._show_called == False: 
            return None
        
        plt.savefig(*self._savefig_args, **self._savefig_kwargs, bbox_inches='tight')
        




    

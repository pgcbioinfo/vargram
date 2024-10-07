"""Main module to generate the mutation profile."""

from . import _profile_renderer
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np
import pandas as pd
from decimal import Decimal


def create_default_colors(num_color, color = '#5E5E5E', single = False):
    """Creates default stack colors."""
    if single:
        return [color]*num_color
    cmap_name = 'viridis'  
    cmap = plt.get_cmap(cmap_name)
    listed_cmap = cmap(np.linspace(0, 1, num_color))
    return [mc.to_hex(color) for color in listed_cmap]


class Profile():

    def __init__(self, data):
        """Initializes Profile attributes."""
        # Class attributes
        self.data = data
        self.plotted_already = False
        self.verbose = False

        # Plot attributes
        self.fig = plt.figure(figsize=[8, 4.8])
        self.closed = False
        self.shown = False

        # Data attributes
        self.x = 'mutation'
        self.y = '' 
        self.ytype = ''
        self.group = 'gene'
        self.stack = 'batch'
        self.key_called = False
        self.threshold = 50
        self.struct = []

        # Aesthetic attributes
        self.xticks_fontsize = 6
        self.xticks_rotation = 90
        self.yticks_fontsize = 7
        self.ylabel = ''
        self.ylabel_fontsize = 'medium'
        self.key_fontsize = 8
        #self.group_label = []
        self.group_fontsize = 'small'
        self.group_title = 'Gene'
        self.stack_label = []
        self.stack_fontsize = 'medium'
        self.stack_color =  ''
        self.stack_title = 'Batch'

    def process(self, **process_kwargs):
        """Creates data for plotting and structure.
        
        Returns
        -------
        pandas.DataFrame
            The processed DataFrame for plotting.

        """
        # Checking if the default columns are present or if not, required arguments are provided
        required_data_attributes = ['x', 'group', 'stack']
        default_data_attributes = [self.x, self.group, self.stack]
        default_columns_detected = all(col in self.data.columns.tolist() for col in default_data_attributes)
        required_attributes_detected = all(arg in process_kwargs.keys() for arg in required_data_attributes) #if arg != 'threshold' and arg != 'ytype')
        if not default_columns_detected and not required_attributes_detected:
            raise ValueError(f"Expected columns or arguments not detected.")

        for process_key in process_kwargs.keys():
            setattr(self, process_key, process_kwargs[process_key])

        # Pivoting dataframe to get x counts
        data_pivoted = self.data.copy()
        if self.y == '': # Choosing what to base counts on
            values_for_counting = 'values_for_counting'
            data_pivoted[values_for_counting] = 1
        elif self.y != '':
            values_for_counting = self.y
        data_pivoted = pd.pivot_table(data_pivoted, index=[self.group, self.x], columns=self.stack, values=values_for_counting, aggfunc="sum", fill_value=0).reset_index() 
        data_pivoted.rename_axis(None, axis=1, inplace = True)

        # Applying threshold, keeping only x
        self.stack_names = [stack for stack in data_pivoted.columns if stack not in [self.group, self.x]]
        if len(self.stack_label) == 0: # Assigning stack_names as labels
            self.stack_label = self.stack_names
        data_filtered = data_pivoted.copy()
        for stack in self.stack_names:
            data_filtered[stack] = data_pivoted[stack].apply(lambda x: 0 if x < self.threshold else x)
        
        # Determining whether to normalize or not
        if self.ytype == '':
            if len(self.stack_names) > 1:
                self.ytype = 'weights'
            else:
                self.ytype = 'counts'
        if self.ylabel == '':
            self.ylabel = self.ytype.title()
        if self.ytype == 'weights':
            for stack in self.stack_names: # Using Decimal for precision
                data_filtered[stack] = data_filtered[stack].apply(lambda x: Decimal(str(x)))
                stack_sum = data_filtered[stack].sum()
                data_filtered[stack] = data_filtered[stack] * Decimal('100') / Decimal(str(stack_sum))
                data_filtered[stack] = data_filtered[stack].apply(lambda x: x.quantize(Decimal('0.01'))) 
        
        # Summing x counts across all stacks
        data_filtered['sum'] = data_filtered[self.stack_names].sum(axis=1)
        if self.ytype == 'weights': # Converting Decimal objects back to float
            for col in self.stack_names + ['sum']:
                data_filtered[col] = data_filtered[col].apply(float)
        data_filtered = data_filtered[data_filtered['sum'] > 0]
        
        # Adding keys if provided
        if self.key_called:
            data_with_keys = pd.merge(data_filtered, self.key_data, on=[self.group, self.x], how='outer')
            data_with_keys.fillna(0, inplace=True)
            if data_filtered['sum'].dtype == np.int64:
                numerical_columns = list(data_with_keys.columns)[2:]
                for col in numerical_columns:
                    data_with_keys[col] = data_with_keys[col].astype(np.int64)
            data_with_keys.reset_index(drop=True, inplace=True)
            self.data_for_plotting = data_with_keys
        else:
            self.data_for_plotting = data_filtered
        data_for_plotting = self.data_for_plotting.copy()
        
        # Getting data for calculating structure
        self.group_names = data_for_plotting[self.group].unique().tolist()
        #if len(self.group_label) == 0: # Assigning group_names as labels
        #    self.group_label = self.group_names 
        unique_counts = []
        for group in self.group_names:
            groups_df = data_for_plotting[data_for_plotting[self.group] == group]
            unique_counts.append(len(groups_df[self.x].unique()))
        
        # Getting data for struct
        self.data_for_struct = pd.DataFrame({self.group: self.group_names, 'count': unique_counts})
        self.data_for_struct.sort_values(by='count', ascending=False, inplace=True)
        self.data_for_struct.reset_index(drop=True, inplace=True)
        
        # Returning for vg.stat()
        if self.data_for_plotting.columns.size > 0 and self.data_for_plotting.shape[0] == 0:
            raise ValueError("Plot DataFrame has no rows. Lowering threshold might help.")
        self.data_for_plotting.sort_values(by=[self.group, self.x], inplace=True)
        self.data_for_plotting.reset_index(drop=True, inplace=True)
        if self.verbose:
            print('** Processed data for plotting. **')
        return self.data_for_plotting
        
    def key(self, **key_kwargs):
        """Obtain the keys."""
        self.key_called = True        
        self.key_data = key_kwargs['key_data']
        self.key_label = key_kwargs['key_labels']
        self.key_color = key_kwargs['key_colors']
        if self.verbose:
            print('** Processed key for profile. **')

    def aes(self, **aes_kwargs):
        for aes_key in aes_kwargs.keys():
            setattr(self, aes_key, aes_kwargs[aes_key])
        if self.verbose:
            print('** Processed aesthetics for profile. **')
    
    def struct_method(self, **struct_kwargs):
        """Obtain the structure for the profile (i.e. the groups to include)."""
        row_groups = struct_kwargs['struct_key'].split('/')
        self.struct = [''.join(col.split()).split(',') for col in row_groups]
        if self.verbose:
            print('** Processed struct. **')

    def plot(self):
        """Create the figure."""
        if self.verbose:
            print('** Plotting **')
        
        # Getting structure of the mutation profile grid
        if len(self.struct) == 0:
            self.struct = _profile_renderer.build_struct(self.data_for_struct, self.group)  
        
        # Creating profile grids
        grids_and_axes = _profile_renderer.build_profile_grid(self.struct, self.data_for_struct, 
                                                      self.group, self.key_called)
        label_grid =  grids_and_axes[0]
        legend_grid =  grids_and_axes[1]
        group_title_axes = grids_and_axes[2]
        barplot_axes = grids_and_axes[3]
        heatmap_axes = grids_and_axes[4]
        
        # Gathering aesthetic attributes
        if self.stack_color == '':
            self.stack_color = create_default_colors(len(self.stack_label))
        x_aes = [self.xticks_fontsize, self.xticks_rotation]
        y_aes = [self.yticks_fontsize, self.ylabel, self.ylabel_fontsize]
        group_aes = [self.group_fontsize, self.group_title]
        stack_aes = [self.stack_fontsize, self.stack_label, self.stack_color, self.stack_title]
        if self.key_called:
            key_aes = [self.key_fontsize, self.key_label, self.key_color]
        else:
            key_aes = [[], [], []]
        
        # Building barplots
        group_labels = []
        _profile_renderer.build_profile(group_title_axes, barplot_axes, heatmap_axes,
                                self.data_for_plotting, self.struct, self.group,
                                self.x, self.fig, self.key_called, 
                                key_aes, self.stack_names, stack_aes,
                                group_labels, x_aes, y_aes)
        # Creating figure y-axis label
        _profile_renderer.build_yaxis_label(self.ylabel, label_grid)
        # Creating figure legend
        _profile_renderer.build_legend(legend_grid, stack_aes, group_aes, group_labels)
        plt.tight_layout()

    def show(self):
        """Displays the generated figure."""
        if not self.plotted_already:
            self.plot()
            self.plotted_already = True
        if self.closed:
            show_fig = plt.figure()
            show_manager = show_fig.canvas.manager
            show_manager.canvas.figure = self.fig
            self.fig.set_canvas(show_manager.canvas)
        plt.show()
        self.shown = True
        if self.verbose:
            print('** Showed figure. **')

    def save(self, **save_kwargs):
        """Saves the generated figure or data."""
        if not self.plotted_already:
            self.plot()
            self.plotted_already = True
        file_extension = save_kwargs['fname'].lower().split('.')[-1]
        # Saving CSV
        if file_extension == 'csv':
            if 'fname' in save_kwargs:
                save_kwargs['path_or_buf'] = save_kwargs.pop('fname')
            self.data_for_plotting.to_csv(**save_kwargs)
            if self.verbose:
                print('** Saved data **')
        # Saving figure
        else:
            if self.shown:
                save_fig = plt.figure()
                save_manager = save_fig.canvas.manager
                save_manager.canvas.figure = self.fig
                self.fig.set_canvas(save_manager.canvas)
            plt.savefig(**save_kwargs, bbox_inches='tight')
            plt.close()
            self.closed = True
            if self.verbose:
                print('** Saved figure **')
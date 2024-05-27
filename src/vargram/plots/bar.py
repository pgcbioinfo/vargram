import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np
from ..methods import bar_wrap, bar_module, stats_module
import pandas as pd

def create_default_colors(num_color, color = '#5E5E5E', single = False):

    if single:
        return [color]*num_color

    cmap_name = 'viridis'  
    cmap = plt.get_cmap(cmap_name)
    listed_cmap = cmap(np.linspace(0, 1, num_color))
    return [mc.to_hex(color) for color in listed_cmap]

class Bar():

    def __init__(self, data):

        # Class attributes
        self.data = data
        self.plotted_already = False

        # Plot attributes
        self.fig = plt.figure(figsize=[8, 4.8])

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

        self.group_label = []
        self.group_fontsize = 'small'
        self.group_title = 'Gene'

        self.stack_label = []
        self.stack_fontsize = 'medium'
        self.stack_color =  ''
        self.stack_title = 'Batch'

    def process(self, **process_kwargs):
        """Creates data for plotting and structure."""

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
            if len(self.stack_label) > 1:
                self.ytype = 'weights'
                self.ylabel = 'Weights'
            else:
                self.ytype = 'counts'
                self.ylabel = 'Counts'
        
        if self.ytype == 'weights':
            for stack in self.stack_names:
                stack_sum = data_filtered[stack].sum()
                data_filtered[stack] = data_filtered[stack] * 100 / stack_sum
        
        # Summing x counts across all stacks
        data_filtered['sum'] = data_filtered[self.stack_names].sum(axis=1)
        data_filtered = data_filtered[data_filtered['sum'] > 0]

        # Getting data for calculating structure
        self.group_names = data_filtered[self.group].unique().tolist()
        if len(self.group_label) == 0: # Assigning group_names as labels
            self.group_label = self.group_names 
        unique_counts = []
        for group in self.group_names:
            groups_df = data_filtered[data_filtered[self.group] == group]
            unique_counts.append(len(groups_df[self.x].unique()))

        self.data_for_struct = pd.DataFrame({self.group: self.group_names, 'count': unique_counts})
        self.data_for_struct.sort_values(by='count', ascending=False, inplace=True)
        self.data_for_struct.reset_index(drop=True, inplace=True)

        # Adding keys if provided
        if not self.key_called:
            self.data_for_plotting = data_filtered
            return None
        
        data_with_keys = pd.merge(data_filtered, self.key_data, on=[self.group, self.x], how='outer')
        data_with_keys.fillna(0, inplace=True)
        data_with_keys.reset_index(drop=True, inplace=True)
        self.data_for_plotting = data_with_keys
        
    def key(self, **key_kwargs):
        self.key_called = True        
        print('** Processed key for barplot. **')

        self.key_data = key_kwargs['key_data']
        self.key_label = key_kwargs['key_labels']
        self.key_color = key_kwargs['key_colors']

    def aes(self, **aes_kwargs):
        print('** Processed aesthetics for barplot. **')

        for aes_key in aes_kwargs.keys():
            setattr(self, aes_key, aes_kwargs[aes_key])
    
    def struct_method(self, **struct_kwargs):
        print('** Processed struct. **')

        row_groups = struct_kwargs['struct_key'].split('/')
        self.struct = [''.join(col.split()).split(',') for col in row_groups]

    def plot(self):
        print('** Plotting **')

        # Getting structure of the mutation profile grid
        if len(self.struct) == 0:
            self.struct = bar_module.build_struct(self.data_for_struct, self.group)      

        # Creating bar grids
        label_grid, legend_grid, group_title_axes, barplot_axes, heatmap_axes  = bar_wrap.build_bar_grid(self.struct, self.data_for_struct, self.group, self.key_called)
        
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
        bar_wrap.build_profile(group_title_axes, barplot_axes, heatmap_axes,
                                self.data_for_plotting, 
                                self.struct,
                                self.group, self.x,
                                self.fig, 
                                self.key_called, key_aes, 
                                stack_aes,
                                group_labels,
                                x_aes, y_aes)
        
        # Creating figure y-axis label
        bar_wrap.build_yaxis_label(self.ylabel, label_grid)

        # Creating figure legend
        bar_wrap.build_legend(legend_grid, stack_aes, group_aes, group_labels)

        plt.tight_layout()

    def stat(self, **stat_kwargs):
        """Saves or displays the generated dataframe
        """

        print('** Saved or displayed stat **')

        if 'print' in stat_kwargs.keys() and stat_kwargs['print'] == True:
            print(self.data_for_plotting.to_markdown())
        else:
            self.data_for_plotting.to_csv(**stat_kwargs)

    def show(self):
        """Displays the generated figure.
        """
        
        print('** Showed figure. **')

        if not self.plotted_already:
            self.plot()
            self.plotted_already = True
        
        plt.show()

    def savefig(self, **savefig_kwargs):
        """Saves the generated figure.
        """

        print('** Saved figure. **')

        if not self.plotted_already:
            self.plot()
            self.plotted_already = True

        plt.savefig(**savefig_kwargs, bbox_inches='tight')

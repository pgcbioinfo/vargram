import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np

def create_default_colors(num_color):
    cmap_name = 'viridis'  
    cmap = plt.get_cmap(cmap_name)
    listed_cmap = cmap(np.linspace(0, 1, num_color))
    return [mc.to_hex(color) for color in listed_cmap]

class Bar():

    def __init__(self, data):

        # Data
        self.data = data
        num_batch_initial = len(self.data['batch'].unique())

        # Plot attributes
        self.fig = plt.figure()

        # Data attributes
        self.x = 'mutation'
        self.xorder = 'position'
        self.y = '' # Indicating the yvalues are counts or weights of unique mutations
        self.group = 'gene'
        self.stack = 'batch'

        # Aesthetic attributes
        self.xticks.fontsize = 6
        self.xticks.rotation = 90

        self.yticks.fontsize = 7
        self.ylabel = ['Instances' if num_batch_initial > 0 else 'Weights']
        self.ylabel.fontsize = 'medium'

        self.key.fontsize = 8

        self.group.fontsize = 'small'
        self.group.label = 'Gene'

        self.stack.fontsize = 'medium'
        self.stack.color =  create_default_colors(num_batch_initial)
        self.stack.label = 'Batch'

    def process(self, **process_kwargs):
        print('Processed data for barplot.')
        pass

    def aes(self, **aes_kwargs):
        print('Processed aesthetics for barplot.')

    def show(self):
        
        print('Showed figure.')
        #plt.tight_layout()
        #plt.show()

    def savefig(self, **savefig_kwargs):

        print('Saved figure.')
        #plt.tight_layout()
        #plt.savefig(**_savefig_kwargs, bbox_inches='tight_layout')
        pass
import matplotlib.pyplot as plt

class Bar():

    def __init__(self, data):

        # Data
        self.data = data

        # Plot attributes
        self.fig = plt.figure()

        # Data attributes
        self.x = 'mutation'
        self.xorder = 'position'
        self.y = 'Instances'
        self.group = 'gene'
        self.stack = 'batch'

        # Aesthetic attributes
        self.xticks.font
        self.xticks.fontsize
        self.xticks.rotation = 90

        self.yticks.font
        self.yticks.fontsize

        self.group.font
        self.group.fontsize
        self.group.label

        self.stack.font 
        self.stack.fontsize
        self.stack.color
        self.stack.label 

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
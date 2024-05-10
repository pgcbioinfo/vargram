import matplotlib.pyplot as plt

class Bar():

    def __init__(self):
        pass

    def process(self, **process_kwargs):
        print('Processed data for barplot.')
        pass

    def aes(self, **aes_kwargs):
        print('Processed aesthetics for barplot.')
        pass

    def show(self):
        
        print('Showed figure.')
        #plt.tight_layout()
        #plt.show()
        pass

    def savefig(self, **savefig_kwargs):

        print('Saved figure.')
        #plt.tight_layout()
        #plt.savefig(**_savefig_kwargs, bbox_inches='tight_layout')
        pass
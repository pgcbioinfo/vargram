"""Module to create test data for mutation profile plot."""

import random
import pandas as pd
import re
from decimal import Decimal

class MyProfileData:

    def __init__(self, key_called, num, ytype):
        self.key_called = key_called
        self.num = num
        self.ytype = ytype

    def create_output(self):
        """Create expected data output."""
        self.output = pd.DataFrame(
            {
                "gene":self.generate_genes(self.num),
                "mutation":self.generate_mutations(self.num)
            }
        )

        # Adding batch counts
        self.num_batches = 2
        self.key_only_ind = random.sample(range(0,  self.num), round(self.num/3))
        for i in range(self.num_batches):
            batch_counts = self.generate_batch_counts(self.num, threshold=50,max_counts=100)
            if self.key_called: 
                # Zero counts in all batches for certain mutations that are only present for the key lineages
                batch_counts = [count if ind not in self.key_only_ind else 0 for ind, count in enumerate(batch_counts)]                
            self.output[f'batch_{i+1}'] = batch_counts.copy()

        # Adding total counts per gene-mutation pair
        self.output['sum'] = self.output.iloc[:, 2:].sum(axis=1)

        # Adding keys
        if  self.key_called:
            self._add_keys_to_output()
        else:
            self.output = self.output[self.output['sum'] != 0]
        
        # Adding position and type
        pos_index = self.output.columns.get_loc('mutation') + 1
        mutation_positions = self.output['mutation'].apply(self.get_position)
        self.output.insert(pos_index, 'position', mutation_positions)
        mutation_types = self.output['mutation'].apply(self.get_type)
        self.output.insert(pos_index+1, 'type', mutation_types)

        # Reordering
        self.output.sort_values(by=['gene', 'position'], inplace=True)
        self.output.reset_index(drop=True, inplace=True)
    
        if self.ytype == 'weights':
            return self.normalized_output()
        else:
            return self.output
    
    def _add_keys_to_output(self):
        """Add key lineages to expected output."""
        # Adding keys
        num_keys = 2

        # Dividing the pure key mutations among the keys
        floor = len(self.key_only_ind) // num_keys
        remainder = len(self.key_only_ind) % num_keys
        num_pure_key_mutations = [floor]*num_keys
        num_pure_key_mutations[-1] = num_pure_key_mutations[-1] + remainder
        for i in range(num_keys):
            # Setting pure key mutations to 0 for this key (which makes it 1 for other keys)
            this_key_ind = self.key_only_ind[:num_pure_key_mutations[i]]
            self.output[f'key_{i+1}'] = [1 if ind not in this_key_ind else 0 for ind in range(self.num)]
            # Removing sampled pure key indices
            self.key_only_ind = self.key_only_ind[num_pure_key_mutations[i]:]
    
    def normalized_output(self):
        """Normalize the mutation counts."""
        self.normalized = self.output.copy()
        numerical_columns = self.normalized.columns.tolist()[4:]
        batches = numerical_columns[:self.num_batches]
        for batch in batches: # Using Decimal for precision
            self.normalized[batch] = self.normalized[batch].apply(lambda x: Decimal(str(x)))
            batch_sum = self.normalized[batch].sum()
            self.normalized[batch] = self.normalized[batch] * Decimal('100') / Decimal(str(batch_sum))
            self.normalized[batch] = self.normalized[batch].apply(lambda x: x.quantize(Decimal('0.01')))        
        self.normalized['sum'] = self.normalized[batches].sum(axis=1)
        for col in numerical_columns:
            self.normalized[col] = self.normalized[col].apply(float)
        return self.normalized

    def create_input(self):
        """Create input that will produce expected output."""
        # Getting batch names
        column_names = list(self.output.columns)
        sum_index = column_names.index("sum")
        type_index = column_names.index("type")
        batch_names = column_names[type_index+1:sum_index]         
        rows = []
        for _, row in self.output.iterrows():
            gene = row['gene']
            mutation = row['mutation'] 
            position = row['position']
            type = row['type']
            # Iterating through batch columns
            for batch in batch_names:
                count = row[batch] 
                # Appending the gene and mutation pair to the rows list 'count' times
                for _ in range(count):
                    rows.append([batch, gene, mutation, position, type])
        # Creating new DataFrame from the rows
        self.input = pd.DataFrame(rows, columns=['batch', 'gene', 'mutation', 'position', 'type'])
        return self.input
    
    def create_keys(self):
        """Create input key lineages."""
        if not self.key_called:
            raise NotImplementedError("Keys cannot be created if key_called is false.")        
        # Getting key names
        column_names = list(self.output.columns)
        sum_index = column_names.index("sum")
        key_names = column_names[sum_index+1:]
        keys = []
        for key_name in key_names:
            keys.append(self.output[self.output[key_name] == 1])
        return keys
    
    def generate_mutations(self, num_mutations, 
                           amino_acids='ACDEFGHIKLMNPQRSTVWY', 
                           max_position=10000):
        """Generate unique mutations."""
        unique_positions = random.sample(range(1, max_position + 1), num_mutations)
        # All mutations are guaranteed to be unique
        mutations = []
        for pos in unique_positions:
            mutation_type = random.randint(0,3)
            if mutation_type <= 1: # substitution
                mutations.append(f'{random.choice(amino_acids)}{pos}{random.choice(amino_acids)}')
            elif mutation_type == 2: # deletion
                mutations.append(f'{random.choice(amino_acids)}{pos}-')
            elif mutation_type == 3: # insertion
                insertion_length = random.randint(1, 10)
                mutations.append(f"{pos}:{''.join(random.choices(amino_acids, k = insertion_length))}") 
        return mutations

    def generate_genes(self, num_genes, num_unique_genes = 10):      
        """Generate genes."""  
        unique_genes = [f'GENE_{i:02}' for i in range(1, num_unique_genes + 1)]
        return [random.choice(unique_genes) for _ in range(num_genes)]

    def generate_batch_counts(self, num_counts, max_counts = 100, threshold = 50):
        """Generate batch counts."""
        counts = [random.randint(threshold, max_counts) for _ in range(num_counts)]
        counts = [count if random.randint(0,1) == 1 else 0 for count in counts]
        return counts

    def get_position(self, mutation):
        """Get position from mutation string."""
        match = re.search(r'(\d+)', mutation)
        return int(match.group())
    
    def get_type(self, mutation):
        """Get type of mutation."""
        if '-' in mutation:
            return 'del'
        elif ':' in mutation:
            return 'in'
        else:
            return 'sub'
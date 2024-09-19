"""Tests whether processed data provided for plotting is correct."""

from vargram import vargram
import random
import pandas as pd
import re
import tempfile
import os
import shutil
import pytest

class MyBarData:

    def __init__(self, key_called, num):

        self.key_called = key_called
        self.num = num

    def create_output(self):
        self.output = pd.DataFrame(
            {
                "gene":self.generate_genes(self.num),
                "mutation":self.generate_mutations(self.num)
            }
        )

        # Adding batch counts
        num_batches = 2
        self.key_only_ind = random.sample(range(0,  self.num), round(self.num/3))
        for i in range(num_batches):
            batch_counts = self.generate_batch_counts(self.num, threshold=50,max_counts=100)
            if self.key_called: # Zero counts in all batches for certain mutations that are only present for the key lineages
                batch_counts = [count if ind not in self.key_only_ind else 0 for ind, count in enumerate(batch_counts)]                
            self.output[f'batch_{i+1}'] = batch_counts.copy()

        # Adding total counts per gene-mutation pair
        self.output['sum'] = self.output.iloc[:, 2:].sum(axis=1)
        if  self.key_called:
            self._add_keys_to_output()
        else:
            self.output = self.output[self.output['sum'] != 0]

        # Reordering
        self.output.sort_values(by=['gene', 'mutation'], inplace=True)
        self.output.reset_index(drop=True, inplace=True)

        return self.output
    
    def _add_keys_to_output(self):

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

    def create_input(self):

        # Getting batch names
        column_names = list(self.output.columns)
        sum_index = column_names.index("sum")
        batch_names = column_names[2:sum_index]            

        rows = []
        for _, row in self.output.iterrows():
            gene = row['gene']
            mutation = row['mutation']
            
            # Iterating through batch columns
            for batch in batch_names:
                count = row[batch]
                
                # Appending the gene and mutation pair to the rows list 'count' times
                for _ in range(count):
                    rows.append([batch, gene, mutation])

        # Creating new DataFrame from the rows
        self.input = pd.DataFrame(rows, columns=['batch', 'gene', 'mutation'])
        mutation_positions = self.input['mutation'].apply(self.extract_position)
        self.input.insert(loc=2, column='position', value=mutation_positions)

        self.input.sort_values(by='position', inplace=True)

        return self.input
    
    def create_keys(self):

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
    
    def generate_mutations(self, num_mutations, amino_acids='ACDEFGHIKLMNPQRSTVWY', max_position=10000):
        unique_positions = random.sample(range(1, max_position + 1), num_mutations)
        sorted_positions = sorted(unique_positions.copy())
        # All mutations are guaranteed to be unique

        mutations = []
        for pos in unique_positions:
            mutation_type = random.randint(0,3)
            if mutation_type <= 1: # substitution
                mutations.append(f'{random.choice(amino_acids)}{pos}{random.choice(amino_acids)}')
            elif mutation_type == 2: # deletion
                pos_ind = sorted_positions.index(pos)
                if pos == max(sorted_positions): 
                    end_position = random.randint(1,10)
                else: # ensuring end_position does not include any unique_position
                    end_positions = [epos for epos in range(sorted_positions[pos_ind+1] - pos + 1) if not any(pos < other_pos < pos + epos for other_pos in sorted_positions)]
                    end_position = pos+random.sample(end_positions, 1)[0]
                mutations.append(f'{pos}-{end_position}')
            elif mutation_type == 3: # insertion
                insertion_length = random.randint(1, 10)
                mutations.append(f'{pos}:{''.join(random.choices(amino_acids, k = insertion_length))}')
        
        return mutations

    def generate_genes(self, num_genes, num_unique_genes = 10):
        
        unique_genes = [f'GENE_{i:02}' for i in range(1, num_unique_genes + 1)]

        return [random.choice(unique_genes) for _ in range(num_genes)]

    def generate_batch_counts(self, num_counts, max_counts = 100, threshold = 50):
        counts = [random.randint(threshold, max_counts) for _ in range(num_counts)]
        counts = [count if random.randint(0,1) == 1 else 0 for count in counts]

        return counts

    def extract_position(self, mutation):

        match = re.search(r'\d+\.?\d*', mutation)
        return float(match.group())
        
@pytest.fixture(params=[(False, 0), (False, 10),(True, 0), (True, 10)])
def bar_data(request):
    key_called, threshold = request.param

    num = random.randint(30,100)
    mbd = MyBarData(key_called=key_called, num=num)
    output = mbd.create_output()
    input = mbd.create_input()

    vg = vargram(data=input)
    vg.bar(threshold=threshold,ytype='counts')

    if key_called == True:
        keys = mbd.create_keys()
        try:
            vargram_test_dir = tempfile.mkdtemp(prefix="vargram_test_dir")
            for i, key in enumerate(keys):
                key_path = os.path.join(vargram_test_dir, f'key_{i+1}.csv')
                key.to_csv(key_path, index=False)
                key.head(5)
                vg.key(key_path)
        finally:
            shutil.rmtree(vargram_test_dir)

    return {"vg":vg, "output":output}

class TestBarData:

    def test_returned_bar_data(self, bar_data):
        """Returned bar data should be equal to saved data. 
        """
        
        vg = bar_data["vg"]
        result = vg.stat()
        try:
            vargram_test_dir = tempfile.mkdtemp(prefix="vargram_test_dir")
            saved_file_path = os.path.join(vargram_test_dir, 'saved_bar_data.csv')
            vg.save(saved_file_path, index=False)
            expected = pd.read_csv(saved_file_path)
        
            assert result.equals(expected)
        finally:
            shutil.rmtree(vargram_test_dir)

    def test_batch_key_sums(self, bar_data):
        """The total sum of all batch and key columns should not be zero.
        """
        vg = bar_data["vg"]
        column_sums = vg.stat().sum(numeric_only=True, axis=1)
        zero_column_sums = column_sums.abs() < 1e-10
        result = zero_column_sums.any()
        expected = False

        assert result == expected
    
    def test_number_mutations(self, bar_data):
        """The length of the dataframe should be equal to the number of unique mutations.
        This is true for the test input created.
        """
        vg = bar_data["vg"]
        output = vg.stat()
        result = len(output['mutation'])
        expected = len(output['mutation'].unique())

        assert result == expected

    def test_bar_data(self, bar_data):
        """The returned bar data should be equal to expected bar data.
        """
        vg = bar_data["vg"]
        expected = bar_data["output"]
        result = vg.stat()

        assert result.equals(expected) 
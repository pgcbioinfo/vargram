"""Tests whether processed data for plotting is correct."""

from vargram import vargram
import random
import pandas as pd
import re
import tempfile
import os
import shutil

def generate_mutations(num_mutations, amino_acids='ACDEFGHIKLMNPQRSTVWY', max_position=1000):
    unique_positions = random.sample(range(1, max_position + 1), num_mutations)
    
    return [f'{random.choice(amino_acids)}{pos}{random.choice(amino_acids)}' for pos in unique_positions]

def generate_genes(num_genes, num_unique_genes = 10):
    
    unique_genes = [f'GENE_{i:02}' for i in range(1, num_unique_genes + 1)]

    return [random.choice(unique_genes) for _ in range(num_genes)]

def generate_batch_counts(num_counts, max_counts = 100, threshold = 50):
    counts = [random.randint(threshold, max_counts) for _ in range(num_counts)]
    counts = [count if random.randint(0,1) == 1 else 0 for count in counts]

    return counts

def extract_position(mutation):

    match = re.search(r'\d+\.?\d*', mutation)
    return float(match.group())

def create_output(num):
    output_df = pd.DataFrame(
        {
            "gene":generate_genes(num),
            "mutation":generate_mutations(num)
        }
    )

    # Adding batch counts
    num_batches = 2
    for i in range(num_batches):
        output_df[f'batch_{i+1}'] = generate_batch_counts(num, threshold=50,max_counts=100)

    # Adding total counts per gene-mutation pair
    output_df['sum'] = output_df.iloc[:, 2:].sum(axis=1)
    output_df = output_df[output_df['sum'] != 0]

    # Reordering
    output_df.sort_values(by=['gene', 'mutation'], inplace=True)
    output_df.reset_index(drop=True, inplace=True)

    return output_df

def create_input(output):

    # Getting batch names
    column_names = list(output.columns)
    batch_names = column_names[2:-1]

    rows = []
    for _, row in output.iterrows():
        gene = row['gene']
        mutation = row['mutation']
        
        # Iterate through batch columns
        for batch_col in batch_names:
            count = row[batch_col]
            
            # Append the gene and mutation pair to the rows list 'count' times
            for _ in range(count):
                rows.append([batch_col, gene, mutation])

    # Create new DataFrame from the rows
    input_df = pd.DataFrame(rows, columns=['batch', 'gene', 'mutation'])
    mutation_positions = input_df['mutation'].apply(extract_position)
    input_df.insert(loc=2, column='position', value=mutation_positions)

    input_df.sort_values(by='position', inplace=True)
    
    return input_df

class TestPlotData:

    def setup_method(self):
        self.num = random.randint(30,100)
        self.output = create_output(self.num)
        self.input = create_input(self.output)

        self.vg = vargram(data=self.input)
        self.vg.bar(threshold=20,ytype='counts')
        self.returned_plot_data = self.vg.stat()

    def test_returned_plot_data(self):
        
        vargram_test_dir = tempfile.mkdtemp(prefix="vargram_test_dir")
        try:
            saved_file_path = os.path.join(vargram_test_dir, 'saved_plot_data.csv')
            self.vg.save(saved_file_path, index=False)
            saved_plot_data = pd.read_csv(saved_file_path)
        
            assert self.returned_plot_data.equals(saved_plot_data)
        finally:
            shutil.rmtree(vargram_test_dir)

    def test_plot_data(self):

        assert self.returned_plot_data.equals(self.output) 
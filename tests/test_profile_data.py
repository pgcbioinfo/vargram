"""Tests whether processed data provided for plotting is correct."""

from create_profile_data import MyProfileData
from vargram import vargram
import matplotlib.pyplot as plt
import random
import pandas as pd
import tempfile
import os
import shutil
import pytest


@pytest.fixture(params=[(False, 0, 'counts'), (False, 10, 'counts'),
                        (False, 0, 'weights'), (False, 10, 'weights'),
                        (True, 0, 'counts'), (True, 10, 'counts'),
                        (True, 0, 'weights'), (True, 10, 'weights')])
def profile_data(request):
    """Feed input to VARGRAM based on different parameters and get output."""
    key_called, threshold, ytype = request.param
    num = random.randint(30,100)

    # Creating test outputs and inputs
    mbd = MyProfileData(key_called=key_called, num=num, ytype=ytype)
    output = mbd.create_output()
    input = mbd.create_input()

    # Getting vargram output
    vg = vargram(data=input, format='_test')
    vg.profile(threshold=threshold,ytype=ytype)
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
    yield {"vg":vg, "output":output}
    plt.close()


class TestProfileData:

    def test_returned_data(self, profile_data):
        """Returned profile data should be equal to saved data. """ 
        vg = profile_data["vg"]
        result = vg.stat()
        try:
            vargram_test_dir = tempfile.mkdtemp(prefix="vargram_test_dir")
            saved_file_path = os.path.join(vargram_test_dir, 'saved_profile_data.csv')
            vg.save(saved_file_path, index=False)
            expected = pd.read_csv(saved_file_path)        
            assert result.equals(expected)
        finally:
            shutil.rmtree(vargram_test_dir)

    def test_batch_key_sums(self, profile_data):
        """The total sum of all batch and key columns should not be zero."""
        vg = profile_data["vg"]
        column_sums = vg.stat().sum(numeric_only=True, axis=1)
        zero_column_sums = column_sums.abs() < 1e-10
        result = zero_column_sums.any()
        expected = False
        assert result == expected
    
    def test_number_mutations(self, profile_data):
        """The length of the dataframe should be equal to the number of unique mutations.
        This is true for the test input created.
        """
        vg = profile_data["vg"]
        output = vg.stat()
        result = len(output['mutation'])
        expected = len(output['mutation'].unique())
        assert result == expected
    
    def test_order_mutations(self, profile_data):
        """Per gene, the positions should only increase as you go down the rows."""
        vg = profile_data["vg"]
        output = vg.stat()
        genes = output['gene'].unique().tolist()

        orders = []
        for gene in genes:
            gene_output = output[output['gene'] == gene]
            positions = list(gene_output['position'])
            position_increasing = all(positions[i] < positions[i + 1] for i in range(len(positions) - 1))
            orders.append(position_increasing)
        
        result = all(orders)
        expected = True
        assert result == expected

    def test_whole_data(self, profile_data):
        """The returned profile data should be equal to expected profile data."""
        vg = profile_data["vg"]
        expected = profile_data["output"]
        result = vg.stat()
        assert result.equals(expected) 
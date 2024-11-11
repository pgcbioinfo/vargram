"""Tests whether profile plot is correct."""

from vargram.plots._profile_renderer import build_struct, build_profile_grid
import random
import pandas as pd


class TestProfileRenderer:

    def setup_method(self):
        # Creating group_counts dataframe
        counts = [50,48,29,10,1,20,2] # Counts are counts of unique mutations per gene
        self.num_genes = len(counts)
        genes = []
        for i in range(self.num_genes):
            genes.append(f'GENE_{i+1}')

        self.gene_counts = pd.DataFrame({
            'my_genes':genes,
            'count':counts
        })

        self.struct = [['GENE_1'], ['GENE_2', 'GENE_5'], ['GENE_3','GENE_4','GENE_7'], ['GENE_6']]
        self.axes = build_profile_grid(self.struct, self.gene_counts, 'my_genes', 'False')

        data_max = max(counts)
        self.predefined_max = 40
        self.max = max(data_max, self.predefined_max)

    def test_struct(self):
        """Asserts that struct algorithm is followed.
        
        Algorithm: Per row, take largest. If len(largest) > max, max = len(largest).
        If the next number added to current sum > max, exclude. Otherwise, include in row.
        Order matters. (In source, the count column is sorted first.)
        """
        expected = self.struct
        result = build_struct(self.gene_counts, 'my_genes', max_per_row=self.predefined_max)

        assert result == expected

    def test_num_axes(self):
        """Number of axes must be equal to number of genes."""

        num_group_title_axes = len(self.axes[2])
        num_barplot_axes = len(self.axes[3])
        num_heatmap_axes = len(self.axes[4])

        num_axes = [num_group_title_axes, num_barplot_axes, num_heatmap_axes]
        equal_to_num_genes = num_axes.count(self.num_genes) == len(num_axes)

        assert equal_to_num_genes == True

    def test_bar_width_ratios(self):
        """The width ratios of each bar ax must sum to the maximum row length (in terms of number of mutations)."""
        sum_width_ratios = []
        for bar_ax in self.axes[3]:
            bar_grid = bar_ax.get_gridspec()
            width_ratios = bar_grid.get_width_ratios()
            sum_width_ratios.append(sum(width_ratios))

        equal_to_max = sum_width_ratios.count(self.max) == len(sum_width_ratios)

        assert equal_to_max == True
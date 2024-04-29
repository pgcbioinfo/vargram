import matplotlib.pyplot as plt
import matplotlib.gridspec as mg
import matplotlib.patches as mp
import matplotlib.text as mt
from . import stats_module 
from . import bar_module 
import copy

def build_bar_grid(gene_rows, grid_width_counts, key_called):
    """Creates the whole GridSpec objects on which to place the plots.

    Args:
        struct (list): Contains the structure of the mutation profile, i.e. gene rows.
        grid_width_counts (pandas.Dataframe): Contains mutation counts per gene to be used for width ratios of grids
        key_called (bool): Determines whether a key lineage was called or not.
    
    Returns:
        label_grid (mg.GridSpecFromSubplotSpec): Grid for the figure y-axis label.
        legend_grid (mg.GridSpecFromSubplotSpec): Grid for the batch (and gene) legend.
        gene_title_axes (list): List of Axes objects on which to plot the gene title.
        gene_mutation_axes (list): List of Axes objects on which to plot the bars.
        gene_key_axes (list): List of Axes objects on which to plot the heatmap for key mutations.
    """
    nrow = len(gene_rows)

    # Main, outermost grid: 1 col for bar ylabel, 1 col for mutation profile, 1 col for legend
    bar_grid = mg.GridSpec(nrow, 3, width_ratios = [0.5, 21, 1])

    # Creating grid for the label and legend columns
    label_grid = mg.GridSpecFromSubplotSpec(1, 1, bar_grid[:, 0])
    legend_grid = mg.GridSpecFromSubplotSpec(2, 1, bar_grid[:, 2], hspace=0.1)

    # Getting maximum width
    width_max = 0
    all_width_ratios = []
    for (i, gene_row) in enumerate(gene_rows):
        width_ratios = grid_width_counts[grid_width_counts['gene'].isin(gene_row)]
        width_ratios = width_ratios['count'].tolist()
        all_width_ratios.append(width_ratios)
        if sum(width_ratios) > width_max:
            width_max = sum(width_ratios)


    # Creating grids for each row of the mutation profile
    modified_gene_rows = copy.deepcopy(gene_rows)
    gene_row_grids = []
    for (i, gene_row) in enumerate(modified_gene_rows):

        row_length = len(gene_row)
        if sum(all_width_ratios[i]) < width_max:
            row_length += 1
            filler_length = width_max-sum(all_width_ratios[i])
            all_width_ratios[i].append(filler_length)
            modified_gene_rows[i].append('filler')

        if key_called:
            height_ratios=[1,7,1.5]
            gene_row_grid = mg.GridSpecFromSubplotSpec(3, row_length, bar_grid[i, 1], 
                                                       width_ratios=all_width_ratios[i], 
                                                       height_ratios=height_ratios,
                                                       wspace=0.1,
                                                       hspace=0.1)
        else:
            height_ratios=[1,3]
            gene_row_grid = mg.GridSpecFromSubplotSpec(2, row_length, bar_grid[i, 1], 
                                                       width_ratios=all_width_ratios[i],
                                                       height_ratios=height_ratios, 
                                                       wspace=0.1,
                                                       hspace=0.1)

        gene_row_grids.append(gene_row_grid)

    gene_title_axes = []
    gene_mutation_axes = []
    gene_key_axes = []

    for (i, gene_row) in enumerate(modified_gene_rows):
        gene_row_grid = gene_row_grids[i]
        for (j, gene) in enumerate(gene_row):
            if gene == 'filler':
                continue
            gene_title_grid = gene_row_grid[0, j]
            gene_mutation_grid = gene_row_grid[1, j]
            
            # Creating subplot for gene titles
            gene_title_ax = plt.subplot(gene_title_grid)

            # Creating subplot for the gene barplots
            if j == 0:
                gene_mutation_ax = plt.subplot(gene_mutation_grid)
            else:
                first_gene_index = len(gene_mutation_axes) - j
                gene_mutation_ax = plt.subplot(gene_mutation_grid, sharey=gene_mutation_axes[first_gene_index])

            gene_title_axes.append(gene_title_ax)
            gene_mutation_axes.append(gene_mutation_ax)

            # Creating subplot for the gene key mutations
            if key_called:
                gene_key_grid = gene_row_grid[2, j]
                gene_key_ax = plt.subplot(gene_key_grid)
                gene_key_axes.append(gene_key_ax)      
    
    return label_grid, legend_grid, gene_title_axes, gene_mutation_axes, gene_key_axes

def build_mutation_profile(gene_title_axes, gene_mutation_axes, gene_key_axes, summary_counts, gene_rows, fig, key_called, key_labels, key_colors, batches, bar_colors, gene_labels):
    """Generates the full mutation profile including labels on the defined grids.

    Args:
        gene_title_axes (list): List of Axes objects on which to plot the gene title
        gene_mutation_axes (list): List of Axes objects on which to plot the bars
        gene_key_axes (list): List of Axes objects on which to plot the heatmap for key mutations
        summary_counts (pandas.Dataframe): The DataFrame containing summary mutation counts per batch including key mutations.
        gene_rows (list): The structure of the mutation profile
        fig (matplotlib.figure.Figure): The Figure object of the entire VARGRAM bar plot.
        key_called (bool): Determines whether a key lineage was called or not.
        key_labels (list): Names of the key lineages provided.
        batches (list): Names of the batches.
        bar_colors (list): Colors of the batches.
        gene_labels (list): List of gene labels that exceed the subplot box.
    
    Returns:
        texts (list): The List of Axes text objects.
        bounds (list): The Axes bounding box object of the texts.
    """

    # Generating colormaps for each key lineage
    heat_cmaps = [] 
    for key_color in key_colors:
        heat_cmap = bar_module.create_colormap(color=key_color)
        heat_cmaps.append(heat_cmap)

    # Getting maximum stacked bar height across each gene
    max_bar_heights = []
    for gene_row in gene_rows:
        max_bar_height = 0
        for gene in gene_row:
            sum_gene = summary_counts[summary_counts['gene'] == gene]
            if max(sum_gene['sum']) > max_bar_height:
                max_bar_height = max(sum_gene['sum'])
        gene_row_max_heights = [max_bar_height]*len(gene_row)
        max_bar_heights += gene_row_max_heights

    # Flattening mutation profile structure and getting first gene on each row
    flattened_genes = [item for gene_row in gene_rows for item in gene_row]
    first_genes = [gene_row[0] for gene_row in gene_rows]

    # Going over each gene and building the barplot, heatmap and gene title
    for (i, gene) in enumerate(flattened_genes):
        sc_gene = summary_counts[summary_counts['gene'] == gene]

        if gene in first_genes:
            suppress_spline = False
            suppress_label = False
        else:
            suppress_spline = True
            suppress_label = True

        # Adding key mutations for gene
        if key_called:
            ax_heat = gene_key_axes[i]
            heat_cmap = bar_module.create_colormap()
            bar_module.build_gene_heatmap(ax_heat, sc_gene, key_labels, heat_cmaps, suppress_label)
        
        # Adding text for gene
        ax_text = gene_title_axes[i]
        bar_module.build_gene_text(ax_text, gene, fig, gene_labels)

        # Creating unit barplot for gene
        ax_bar = gene_mutation_axes[i]
        floor = [0]*len(sc_gene)
        for (batch, color) in zip(batches, bar_colors):
            bar_module.build_gene_barplot(ax_bar, sc_gene['mutation'], sc_gene[batch], floor=floor,
                                          batch_color=color, suppress_spline=suppress_spline, key_called=key_called, max_height=max_bar_heights[i])
            floor += sc_gene[batch]            

def build_yaxis_label(label, label_grid):
    """Generates the label of the figure y-axis.
    
    Args:
        label (string): The y-axis label.
        label_grid (mg.GridSpecFromSubplotSpec): Grid for the figure y-axis label.

    Returns:
        None.
    """
    # text() settings
    fontsize='medium'
    ax_label = plt.subplot(label_grid[:, 0])

    # Creating label
    xlims = ax_label.get_xlim()
    ylims = ax_label.get_ylim()
    ax_label.text(xlims[1]/2, ylims[1]/2, label, ha='center', va='center', transform=ax_label.transAxes, fontsize=fontsize,rotation=90)

    # Removing spines and ticks
    ax_label.set_yticks([])
    ax_label.set_xticks([])
    ax_label.spines["top"].set_visible(False)
    ax_label.spines["bottom"].set_visible(False)
    ax_label.spines["left"].set_visible(False)
    ax_label.spines["right"].set_visible(False)




def build_legend(legend_grid, batches, batch_colors, legend_title, gene_labels):
    """Generates the label of the figure y-axis.
    
    Args:
        legend_grid (mg.GridSpecFromSubplotSpec): Grid for the batch (and gene) legend.
        batches (list): The list of batches.
        batch_colors (lists): Colors of the batches.
        legend_title (string): Title of the batch legend.
        gene_labels (list): List of gene labels.

    Returns:
        None.
    """

    # legend() settings
    if len(gene_labels) == 0:
        ax_batch_legend = plt.subplot(legend_grid[:, 0])
    else:
        ax_batch_legend = plt.subplot(legend_grid[0, 0])
        ax_gene_legend = plt.subplot(legend_grid[1, 0])
    batch_title=legend_title
    gene_title='Gene'
    fontsize='medium'
    frameon=False
    alignment='left'
    batch_loc='lower left'
    gene_loc='upper left'
    

    # Setting batch legend handles
    batch_legend_handles = [mp.Patch(color=color, label=label) for color, label in zip(batch_colors, batches)]

    # Creating batch legend
    bl = ax_batch_legend.legend(handles=batch_legend_handles, title=batch_title, fontsize=fontsize, frameon=frameon, alignment=alignment, loc=batch_loc, bbox_to_anchor=(-0.5, 0), borderaxespad=0)

    # Removing batch ax spines and ticks
    bar_module.spine_remover(ax_batch_legend)

    if len(gene_labels) != 0:

        # Creating handles
        legend_handles = [Text(str(i+1)) for i in range(len(gene_labels))]

        # Creating gene legend
        gl = ax_gene_legend.legend(legend_handles, gene_labels, title=gene_title, fontsize=fontsize, frameon=frameon, alignment=alignment, loc=gene_loc, bbox_to_anchor=(-0.5, 1), borderaxespad=0,
           handler_map={handle: TextHandler() for handle in legend_handles},
           handletextpad=0.1,
           labelspacing=1)
        
        # Removing gene ax spines and ticks
        bar_module.spine_remover(ax_gene_legend)

class Text(object):
    def __init__(self, text):
        self.text = text

class TextHandler(object):
    def legend_artist(self, legend, text_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = mt.Text(x=width/4, y=0, text=text_handle.text, bbox=dict(facecolor='none', boxstyle='Square'),
                        weight='bold',
                        verticalalignment=u'baseline', 
                        horizontalalignment=u'left', multialignment=None, 
                        fontproperties=None, linespacing=None, 
                        rotation_mode=None)
        handlebox.add_artist(patch)

        return patch
    





    
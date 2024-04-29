"""Module to generate a bar plot of unique mutations and their counts, i.e. the VARGRAM."""

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import pandas as pd
import numpy as np

def build_struct(gene_counts, max_per_row = 40):
    """Determines the optimum structure of the genes in the bar plot.

    Args:
        gene_counts (pandas.DataFrame): The DataFrame containing genes and their unique no. of mutations. 
                                        Possibly includes key mutations.
        max_per_row (int): Initial maximum number of mutations per row.
    
    Returns:
        struct (list): The structure of the plot where each row gives the list of genes for that row.
    """

    gg = gene_counts['gene'].tolist()
    cc = gene_counts['count'].tolist()


    struct = [] # list of genes per row
    struct_len = [] # list of total mutation counts per row
    while len(gg) > 0:
        # Each iteration determines the genes for a row

        # Setting length of largest gene as max_per_row
        # Since it is larger, it takes its own row
        largest_count = max(cc)
        if largest_count >= max_per_row: 
            max_per_row = largest_count

            largest_index = cc.index(largest_count)
            largest_gene = gg[largest_index]
            struct.append([largest_gene])

            cc.remove(largest_count)
            gg.remove(largest_gene)

            continue

        # Of the remaining genes, find all that sum to less than or equal to max_per_row
        gene_row = [gg[0]]
        current_sum = cc[0]
        for (i, gene) in enumerate(gg):
            if i == 0:
                continue
            if current_sum + cc[i] <= max_per_row:
                gene_row += [gene]
                current_sum += cc[i]
        struct.append(gene_row)

        # Remove these genes
        indices_to_remove = [i for i, gene in enumerate(gg) if gene in gene_row]
        gg = [gene for i, gene in enumerate(gg) if i not in indices_to_remove]
        cc = [count for i, count in enumerate(cc) if i not in indices_to_remove]

    
    return struct

def compute_gridsize(struct, nkey):
    """Determines the rows and columns of the GridSpec for the bar plot.

    Args:
        struct (pandas.DataFrame): The DataFrame containing genes and their no. of mutations.
        nkey (integer): The number of reference lineages.
    
    Returns:
        nrow (integer): The number of rows. 
        ncol (integer): The number of columns.
    """
    # Counting columns
    ncol = 3 # one for y-axis label, one for the main plot, one for the legend

    # Setting initial number of rows 
    # gs_row = (2 + No. of Ref. Lineages)*(No. of Gene Rows)
    nrow = (2 + nkey)*len(struct)

    return nrow, ncol

def build_gene_barplot(ax_bar, categories, heights, floor, batch_color, suppress_spline, key_called, max_height):
    """Generates the individual barplot of a gene.

    Args:
        ax_bar (matplotlib.axes.Axes): The subplot where to place the barplot for a particular gene.
        categories (pandas.DataFrame): All mutations for a particular gene.
        heights
        floor (list): The bottom y-values at which to start plotting the bars.
        batch_color (str): The color for this particular batch.
        suppress_spline (bool): Determines whether the y-axis and left spine will be shown.
        max_height (int): Maximum height of a stacked bar.

    Returns:
        None.
    """

    # bar() settings
    width = 0.75
    edgecolor = 'black'
    linewidth = 1

    # Creating barplot
    ax_bar.bar(x = categories,
           height = heights,
           bottom = floor,
           color = batch_color,
           width = width,
           edgecolor = edgecolor,
           linewidth = linewidth
           )
    
    # Removing spines
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.spines["bottom"].set_visible(False)

    # Adjusting limits of x-axis and y-axis
    ax_bar.set_xlim(-0.5, len(categories) - 0.5)
    if max_height != 0: # Avoids UserWarning
        ax_bar.set_ylim(0.0, max_height+1) 
    else:
        ax_bar.set_ylim(0.0, 1.0)

    # Removing x-axis labels if key lineage is called
    if key_called:
        ax_bar.xaxis.set_visible(False)
    else:
        ax_bar.tick_params(axis='x', rotation=90, labelsize=6)

    # Leaving y-axis and left spine depending on whether this corresponds to first gene on the bar row
    if suppress_spline or max_height == 0:
        ax_bar.yaxis.set_visible(False)
        ax_bar.spines['left'].set_visible(False)
    else:
        ax_bar.spines['left'].set_linewidth(1.5)
        ax_bar.spines['left'].set_position(('outward', 5)) 
        ax_bar.yaxis.set_tick_params(width=1.5, labelsize=7)
        yticks = np.linspace(0, max_height, 5)
        yticks = np.round(yticks).astype(int)
        ax_bar.set_yticks(yticks[yticks != 0])    

    return None

def create_colormap(cmap_name = "lineage", color = "#5E5E5E"):
    """Generate the colormap of a reference lineage.
    
    Args:
        cmap_name (str): The name of the reference lineage.
        color (str): Hex code of color to use for this lineage.
    
    Returns:
        matplotlib.colors.LinearSegmentedColormap: The colormap.
    """

    cmap_colors = ["#D5D5D5", color]

    return mc.LinearSegmentedColormap.from_list(cmap_name, cmap_colors)

def build_gene_heatmap(ax_heat, gene_mutations, key_labels, cmaps, suppress_label):
    """Generates the individual heatmap of a reference lineage.

    Args:
        ax_heat (matplotlib.axes.Axes): The subplot where to place the heatmap for a particular gene.
        gene_mutations (pandas.DataFrame): Summary counts for a particular gene including key mutations.
        key_labels (string): The name of the reference lineage.
        cmaps (matplotlib.colors.Colormap): List of colormaps to use.
        suppress_label (bool): Determines whether key lineage label should be shown.

    Returns:
        None. 
    """
    # Converting mutation names into binaries for imshow()
    mutations_matrix = []
    for key_label in key_labels:
        mutations_in_binary = gene_mutations[key_label]
        mutations_matrix.append(mutations_in_binary)

    # imshow() settings
    heatmap_border_color = 'black'
    heatmap_border_linewidth = 3
    heatmap_partition_color = 'white'
    heatmap_partition_linewidth = 1.5

    # Creating heatmap
    reversed_cmaps=cmaps[::-1]
    mutation_names = gene_mutations['mutation']
    for i, row in enumerate(reversed(mutations_matrix)):
        ax_heat.imshow([row], cmap=reversed_cmaps[i], vmin=0, vmax=1, extent=(-0.5, len(mutation_names)-0.5, i-0.5, i+0.5), aspect='auto')
    ax_heat.tick_params(axis='x', rotation=90, labelsize=6)
    ax_heat.set_xticks(np.arange(len(mutation_names)))
    ax_heat.set_xticklabels(mutation_names)
    ax_heat.vlines(x=np.arange(0, len(mutation_names)-1)+0.5, ymin = -0.5, ymax = len(key_labels)-0.5, color = heatmap_partition_color, linewidth = heatmap_partition_linewidth)

    # Adding key lineage label
    if suppress_label:
        ax_heat.set_yticks([])
    else:
        ax_heat.set_yticks(list(range(len(key_labels))))
        ax_heat.set_yticklabels(key_labels[::-1])
        ax_heat.yaxis.set_tick_params(labelsize=8)

    # Creating heatmap border
    ax_heat.vlines(x=-0.5, ymin=-0.5, ymax=len(key_labels)-0.5, color = heatmap_border_color, linewidth = heatmap_border_linewidth)
    ax_heat.vlines(x=len(mutation_names)-0.5, ymin=-0.5, ymax=len(key_labels)-0.5, color = heatmap_border_color, linewidth = heatmap_border_linewidth)
    ax_heat.hlines(y=np.linspace(-0.5, len(key_labels)-0.5, len(key_labels)+1), xmin= -0.5, xmax=len(mutation_names)-0.5, color = heatmap_border_color, linewidth = heatmap_border_linewidth)

    return None

def build_gene_text(ax_text, gene_name, fig_text, gene_labels):
    """Generates the gene label above the barplot.

    Args:
        ax_text (matplotlib.axes.Axes): The subplot where to place the gene name text for a particular gene.
        gene_name (str): The text.
        fig_text (matplotlib.figure.Figure): The Figure object of the entire VARGRAM bar plot.
        gene_labels (list): List of gene labels that exceed the subplot box.

    Returns:
        text ().
    """

    # text() settings
    fontsize = 'small'
    weight = 'bold'

    # Creating text
    xlims = ax_text.get_xlim()
    ylims = ax_text.get_ylim()
    t = ax_text.text(xlims[1]/2, ylims[1]/2, gene_name, ha='center', va='center', transform=ax_text.transAxes, fontsize=fontsize, weight=weight)

    # Removing spines and ticks
    ax_text.set_yticks([])
    ax_text.set_xticks([])
    ax_text.spines['left'].set_linewidth(1.5)
    ax_text.spines["top"].set_linewidth(1.5)
    ax_text.spines["right"].set_linewidth(1.5)
    ax_text.spines["bottom"].set_linewidth(1.5)

    # Determining if text exceeds its subplot
    b = t.get_window_extent(renderer=fig_text.canvas.get_renderer()).transformed(ax_text.transData.inverted())

    if b.height > 0.8:
        margin = 0.5

        # Get the current figsize
        fig_width, fig_height = ax_text.figure.get_size_inches()

        # Update the figsize to accommodate the required margin
        new_fig_height = fig_height + margin
        ax_text.figure.set_size_inches((4/3)*new_fig_height, new_fig_height, forward=True)

    if b.width > 1.2:
        margin = 0.75
        width_margin = margin

        # Get the current figsize
        fig_width, fig_height = ax_text.figure.get_size_inches()

        # Update the figsize to accommodate the required margin
        new_fig_width = fig_width + width_margin
        ax_text.figure.set_size_inches(new_fig_width, (3/4)*new_fig_width, forward=True)
    
    # Rechecking
    b = t.get_window_extent(renderer=fig_text.canvas.get_renderer()).transformed(ax_text.transData.inverted())
    if b.width > 1.2:
        no_gene_labels = len(gene_labels)
        gene_labels.append(gene_name)
        t.set_text('{}'.format(no_gene_labels+1))


def spine_remover(ax):
    """Removes the spine of an Axes object .
    
    Args:
        ax (matplotlib.Axes): The Axes.

    Returns:
        None.
    """

    ax.set_yticks([])
    ax.set_xticks([])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)

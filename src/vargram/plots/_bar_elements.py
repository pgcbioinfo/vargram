"""Module to generate a bar plot of unique mutations and their counts, i.e. the VARGRAM."""

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import pandas as pd
import numpy as np

def check_xticks_overlap(xticks):
    """Checks whether there is an overlap between the x-axis tick labels.

    Args:
        xticks (matplotlib.axes.Axes.get_xticklabels): X-axis tick labels.

    Returns: 
        bool: True if there is an overlap, False otherwise.
    """

    bb1 = xticks[0].get_window_extent()
    bb2 = xticks[1].get_window_extent()

    return not (bb1.x1 < bb2.x0 or bb2.x1 < bb1.x0 or bb1.y1 < bb2.y0 or bb2.y1 < bb1.y0)


def build_group_barplot(ax_bar, categories, heights, floor, batch_color, suppress_spline, key_called, max_height, x_aes, y_aes):
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
    x_fontsize = x_aes[0]
    y_fontsize = y_aes[0]

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
        ax_bar.tick_params(axis='x', rotation=90, labelsize=x_fontsize)
        
        # Checking for x-axis overlap
        xticks = ax_bar.get_xticklabels()
        if len(xticks) > 1:
            margin = 0.1
            while check_xticks_overlap(xticks):

                # Get the current figsize
                fig_width, fig_height = ax_bar.figure.get_size_inches()

                # Update the figsize to accommodate the required margin
                new_fig_width = fig_width + margin
                ax_bar.figure.set_size_inches((new_fig_width,(3/4)*new_fig_width), forward=False)
            

    # Leaving y-axis and left spine depending on whether this corresponds to first gene on the bar row
    if suppress_spline or max_height == 0:
        ax_bar.yaxis.set_visible(False)
        ax_bar.spines['left'].set_visible(False)
    else:
        ax_bar.spines['left'].set_linewidth(1.5)
        ax_bar.spines['left'].set_position(('outward', 5)) 
        ax_bar.yaxis.set_tick_params(width=1.5, labelsize=y_fontsize)
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

    # Checking for x-axis overlap
    xticks = ax_heat.get_xticklabels()
    if len(xticks) > 1:
        margin = 0.1
        while check_xticks_overlap(xticks):

            # Get the current figsize
            fig_width, fig_height = ax_heat.figure.get_size_inches()

            # Update the figsize to accommodate the required margin
            new_fig_width = fig_width + margin
            ax_heat.figure.set_size_inches((new_fig_width,(3/4)*new_fig_width), forward=False)

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

def build_group_text(ax_text, gene_name, fig_text, gene_labels):
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
    b = t.get_window_extent(renderer=ax_text.figure.canvas.get_renderer()).transformed(ax_text.transData.inverted())

    if b.height > 0.8:
        margin = 0.5

        # Get the current figsize
        fig_width, fig_height = ax_text.figure.get_size_inches()

        # Update the figsize to accommodate the required margin
        new_fig_height = fig_height + margin
        ax_text.figure.set_size_inches(((4/3)*new_fig_height, new_fig_height), forward=False)

    if b.width > 1.2:
        margin = 0.75
        width_margin = margin

        # Get the current figsize
        fig_width, fig_height = ax_text.figure.get_size_inches()

        # Update the figsize to accommodate the required margin
        new_fig_width = fig_width + width_margin
        ax_text.figure.set_size_inches((new_fig_width, (3/4)*new_fig_width), forward=False)
    
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
"""Module that generates the bar structure, grid and plots its elements."""

from . import _profile_elements
import matplotlib.pyplot as plt
import matplotlib.gridspec as mg
import matplotlib.patches as mp
import matplotlib.text as mt
import copy


def build_ordered_struct(group_counts, group_attr, ordered_genes, flat=False, max_per_row=40):
    """Builds the structure based on the genes' start positions

    Parameters
    ----------
    group_counts : pandas.DataFrame
        The DataFrame containing groups and their unique no. of x data.
    group_attr : str
        The group data attribute.
    ordered_genes: list
        A list of ordered genes based on CDS start position.
    flat : bool, default : False
        Determines if groups should be plotted on only one row ("flat") or not ("compact").
    
    Returns
    -------
    list
        The structure of the plot where each row gives the list of groups for that row.
    """
    if flat is True:
        struct = [[gene for gene in ordered_genes]]
        return struct
    
    gg = group_counts[group_attr].tolist() # Gene names
    cc = group_counts['count'].tolist() # Gene counts
    ref_row_length = max(max_per_row, max(cc))
    struct = []
    row = []
    num_row = 0
    row_sum = 0
    for (i, gene) in enumerate(ordered_genes):
        gene_index = gg.index(gene)
        gene_count = cc[gene_index]

        if i == 0: # First gene is the first value in struct, no computation needed
            row.append(gene)
            row_sum += gene_count
            continue

        if row_sum < ref_row_length and num_row == 0:
            row.append(gene)
            row_sum += gene_count
        elif row_sum + gene_count <= ref_row_length:
            row.append(gene)
            row_sum += gene_count
        else:
            num_row += 1
            if row_sum > ref_row_length:
                ref_row_length = row_sum
            struct.append(row)
            row_sum = gene_count
            row = []
            row.append(gene)
    struct.append(row)

    return struct


def build_struct(group_counts, group_attr, flat=False, max_per_row=40):
    """Determines the optimum structure of the groups in the bar plot.

    Parameters
    ----------
    group_counts : pandas.DataFrame
        The DataFrame containing groups and their unique no. of x data.
    group_attr : str
        The group data attribute.
    flat : bool, default False
        Determines if groups should be plotted on only one row (if True).
    max_per_row : int, default 40
        Initial maximum number of x data per row.
    
    Returns
    -------
    list
        The structure of the plot where each row gives the list of groups for that row.

    """
    gg = group_counts[group_attr].tolist()
    cc = group_counts['count'].tolist()
    struct = [] # list of groups per row

    if flat:
        paired_counts = list(zip(gg, cc))
        descending_paired = sorted(paired_counts, key=lambda x: x[1], reverse=True)
        struct = [[group for group, _ in descending_paired]]
        return struct

    while len(gg) > 0: # Each iteration determines the groups for a row
        # Setting length of largest group as max_per_row
        # Since it is larger, it takes its own row
        largest_count = max(cc)
        if largest_count >= max_per_row: 
            max_per_row = largest_count
            largest_index = cc.index(largest_count)
            largest_group = gg[largest_index]
            struct.append([largest_group])
            cc.remove(largest_count)
            gg.remove(largest_group)
            continue

        # Of the remaining groups, find all that sum to less than or equal to max_per_row
        group_row = [gg[0]]
        current_sum = cc[0]
        for (i, group) in enumerate(gg):
            if i == 0:
                continue
            if current_sum + cc[i] <= max_per_row:
                group_row += [group]
                current_sum += cc[i]
        struct.append(group_row)

        # Remove these groups
        indices_to_remove = [i for i, group in enumerate(gg) if group in group_row]
        gg = [group for i, group in enumerate(gg) if i not in indices_to_remove]
        cc = [count for i, count in enumerate(cc) if i not in indices_to_remove]

    return struct


def build_profile_grid(struct, grid_width_counts, group_attr, key_called):
    """Creates the whole GridSpec objects on which to place the plots.

    Parameters
    ----------
    struct : list
        Contains the structure of the profile, i.e. group rows.
    grid_width_counts : pandas.Dataframe
        Contains x counts per group to be used for width ratios of grids
    group_attr : str
        The group attribute of the data
    key_called : bool
        Determines whether a key lineage was called or not.
    
    Returns
    -------
    mg.GridSpecFromSubplotSpec
        Grid for the figure y-axis label.
    mg.GridSpecFromSubplotSpec
        Grid for the stack (and group) legend.
    list
        List of Axes objects on which to plot the group title.
    list
        List of Axes objects on which to plot the bars.
    list
        List of Axes objects on which to plot the heatmap for key data.

    """
    # Main, outermost grid: 1 col for bar ylabel, 1 col for profile, 1 col for legend
    nrow = len(struct)
    bar_grid = mg.GridSpec(nrow, 3, width_ratios=[0.15, 21, 0.5])

    # Creating grid for the label and legend columns
    label_grid = mg.GridSpecFromSubplotSpec(1, 1, bar_grid[:, 0])
    legend_grid = mg.GridSpecFromSubplotSpec(2, 1, bar_grid[:, 2], hspace=0.1)

    # Getting maximum width
    width_max = 0
    all_width_ratios = []
    for (i, group_row) in enumerate(struct):
        group_row_counts = grid_width_counts[grid_width_counts[group_attr].isin(group_row)]
        ordered_group_row_counts = group_row_counts.set_index(group_attr).reindex(group_row).reset_index()
        width_ratios = ordered_group_row_counts['count'].tolist()
        all_width_ratios.append(width_ratios)
        if sum(width_ratios) > width_max:
            width_max = sum(width_ratios)

    # Creating grids for each row of the profile
    modified_group_rows = copy.deepcopy(struct)
    group_row_grids = []
    for (i, group_row) in enumerate(modified_group_rows):
        # Determining if row is not full
        # Then, add a filler group
        row_length = len(group_row)
        if sum(all_width_ratios[i]) < width_max:
            row_length += 1
            filler_length = width_max-sum(all_width_ratios[i])
            all_width_ratios[i].append(filler_length)
            modified_group_rows[i].append('filler')
        
        # Creating row grids
        if key_called:
            height_ratios = [1,7,1.5]
            group_row_grid = mg.GridSpecFromSubplotSpec(3, row_length, bar_grid[i, 1], 
                                                       width_ratios=all_width_ratios[i], 
                                                       height_ratios=height_ratios,
                                                       wspace=0.1,
                                                       hspace=0.1)
        else:
            height_ratios = [1,8.5]
            group_row_grid = mg.GridSpecFromSubplotSpec(2, row_length, bar_grid[i, 1], 
                                                       width_ratios=all_width_ratios[i],
                                                       height_ratios=height_ratios, 
                                                       wspace=0.1,
                                                       hspace=0.1)
        group_row_grids.append(group_row_grid)

    # Creating the axes
    group_title_axes = []
    group_x_axes = []
    group_key_axes = []
    for (i, group_row) in enumerate(modified_group_rows):
        group_row_grid = group_row_grids[i]
        for (j, group) in enumerate(group_row):
            if group == 'filler':
                continue
            group_title_grid = group_row_grid[0, j]
            group_x_grid = group_row_grid[1, j]
            
            # Creating subplot for group titles
            group_title_ax = plt.subplot(group_title_grid)

            # Creating subplot for the group barplots
            if j == 0:
                group_x_ax = plt.subplot(group_x_grid)
            else:
                first_group_index = len(group_x_axes) - j
                group_x_ax = plt.subplot(group_x_grid, sharey=group_x_axes[first_group_index])

            group_title_axes.append(group_title_ax)
            group_x_axes.append(group_x_ax)

            # Creating subplot for the group key x data
            if key_called:
                group_key_grid = group_row_grid[2, j]
                group_key_ax = plt.subplot(group_key_grid)
                group_key_axes.append(group_key_ax)      
    
    return label_grid, legend_grid, group_title_axes, group_x_axes, group_key_axes


def build_profile(group_title_axes, group_bar_axes, group_key_axes, 
                  barplot_data, struct, group_attr, x_attr, fig, aspect, 
                  key_called, key_aes, stack_names, stack_aes, 
                  group_labels, x_aes, y_aes):
    """Generates the full profile including labels on the defined grids.

    Parameters
    ----------
    group_title_axes : list
        List of Axes objects on which to plot the group title.
    group_bar_axes : list
        List of Axes objects on which to plot the bars.
    group_key_axes : list
        List of Axes objects on which to plot the heatmap for key data.
    barplot_data : pandas.Dataframe
        The DataFrame containing summary x counts per stack including key values.
    struct : list
        The structure of the barplot/mutation profile.
    group_attr : str
        The column of the groups.
    x_attr : str
        The column of the x values.
    fig : matplotlib.figure.Figure
        The Figure object of the entire VARGRAM bar plot.
    aspect : float
        The aspect ratio (width / height) of the Figure.
    key_called : bool
        Determines whether a key was called or not.
    key_aes : list
        Aesthetic attributes of the key heatmap.
    stack_names : list
        The names of the stacks (from the data provided).
    stack_aes : list
        Aesthetic attributes of the stacks including 
        stack labels (may be user-provided) to put on the legend.
    group_labels : list
        List of group labels (may be user-provided) that exceed the subplot box.
    x_aes : list
        Aesthetic attributes of the x ticks and labels
    y_aes : list
        Aesthetic attributes of the y label and the y label itself
    
    Returns
    -------
    None

    """
    # Defining aesthetic attributes
    stacks = stack_names
    stack_colors = stack_aes[2]
    key_fontsize = key_aes[0]
    key_labels = key_aes[1]
    key_colors = key_aes[2]

    # Generating colormaps for each key lineage
    heat_cmaps = [] 
    for key_color in key_colors:
        heat_cmap = _profile_elements.create_colormap(color=key_color)
        heat_cmaps.append(heat_cmap)

    # Getting maximum stacked bar height across each group
    max_bar_heights = []
    for group_row in struct:
        max_bar_height = 0
        for group in group_row:
            group_barplot_data = barplot_data[barplot_data[group_attr] == group]
            if max(group_barplot_data['sum']) > max_bar_height:
                max_bar_height = max(group_barplot_data['sum'])
        group_row_max_heights = [max_bar_height]*len(group_row)
        max_bar_heights += group_row_max_heights

    # Flattening structure and getting first group on each row
    flattened_groups = [item for group_row in struct for item in group_row]
    first_groups = [group_row[0] for group_row in struct]

    # Going over each group and building the barplot, heatmap and group title
    for (i, group) in enumerate(flattened_groups):
        group_barplot_data = barplot_data[barplot_data[group_attr] == group]

        # Determining whether to suppress splines
        if group in first_groups:
            suppress_spline = False
            suppress_label = False
        else:
            suppress_spline = True
            suppress_label = True

        # Adding key x data for group
        if key_called:
            ax_heat = group_key_axes[i]
            heat_cmap = _profile_elements.create_colormap()
            _profile_elements.build_group_heatmap(ax_heat, group_barplot_data, key_labels, key_fontsize, heat_cmaps, suppress_label)
        
        # Adding text for group
        ax_text = group_title_axes[i]
        _profile_elements.build_group_text(ax_text, group, fig, aspect, group_labels)

        # Creating unit barplot for group
        ax_bar = group_bar_axes[i]
        floor = [0]*len(group_barplot_data)
        for (batch, color) in zip(stacks, stack_colors):
            _profile_elements.build_group_barplot(ax_bar, group_barplot_data[x_attr], group_barplot_data[batch], floor,
                                          color, suppress_spline, key_called, max_bar_heights[i],
                                          x_aes, y_aes)
            floor += group_barplot_data[batch]            


def build_yaxis_label(label, label_grid):
    """Generates the label of the figure y-axis.
    
    Parameters
    ----------
    label : str
        The y-axis label.
    label_grid : mg.GridSpecFromSubplotSpec
        Grid for the figure y-axis label.
    
    Returns
    -------
    None

    """
    # text() settings
    fontsize='large'
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


def build_legend(legend_grid, stack_aes, group_aes, group_labels):
    """Generates the label of the figure y-axis.
    
    Parameters
    ----------
    legend_grid : mg.GridSpecFromSubplotSpec
        Grid for the stack (and group) legend.
    stack_aes : list
        Aesthetic attributes of the stack legend.
    group_aes : str
        Aesthetic attributes of the group legend.
    group_labels : list
        List of group labels.

    Returns
    -------
    None

    """
    # legend() settings
    if len(group_labels) == 0:
        ax_batch_legend = plt.subplot(legend_grid[:, 0])
    else:
        ax_batch_legend = plt.subplot(legend_grid[0, 0])
        ax_group_legend = plt.subplot(legend_grid[1, 0])
    stack_fontsize = stack_aes[0]
    stack_label = stack_aes[1]
    stack_color = stack_aes[2]
    stack_title = stack_aes[3]
    group_fontsize = group_aes[0]
    group_title = group_aes[1]
    frameon=False
    alignment='left'
    if len(group_labels) == 0:
        stack_loc='center'
        bbox_anchor = (0.5, 0.5)
    else:
        stack_loc='lower left'
        bbox_anchor = (-0.25,0)#(-0.5, 0)
    group_loc='upper left'

    # Setting batch legend handles
    batch_legend_handles = [mp.Patch(color=color, label=label) for color, label in zip(stack_color, stack_label)]

    # Creating batch legend
    ax_batch_legend.legend(handles=batch_legend_handles, 
                           title=stack_title, 
                           title_fontsize=stack_fontsize,
                           fontsize=stack_fontsize, 
                           frameon=frameon, 
                           alignment=alignment, 
                           loc=stack_loc, 
                           bbox_to_anchor=bbox_anchor, 
                           borderaxespad=0)

    # Removing batch ax spines and ticks
    _profile_elements.spine_remover(ax_batch_legend)

    if len(group_labels) != 0:
        # Creating handles
        legend_handles = [Text(str(i+1)) for i in range(len(group_labels))]

        # Creating group legend
        ax_group_legend.legend(legend_handles, 
                               group_labels, 
                               title=group_title,
                               title_fontsize=stack_fontsize, 
                               fontsize=group_fontsize, 
                               frameon=frameon, 
                               alignment=alignment, 
                               loc=group_loc, 
                               bbox_to_anchor=(-0.25,1),#(-0.5, 1), 
                               borderaxespad=0,
                               handler_map={handle: TextHandler() for handle in legend_handles},
                               handletextpad=0.5,
                               labelspacing=0.9)
        
        # Removing group ax spines and ticks
        _profile_elements.spine_remover(ax_group_legend)


class Text(object):

    def __init__(self, text):
        self.text = text


class TextHandler(object):
    """Creates the text patch for the group legend."""

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
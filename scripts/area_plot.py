from warnings import Warning

from matplotlib import rcParams
from matplitlib import cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

"""
A set of python scripts to generate an area plot.
"""
# Sets up the matplotlib paramters so that we can save to be edited in 
# illustator if a direct conversion is required. Because it just makes life
# better
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Sets up the order of colors to be used in joint plots.
colors_order = ['Reds', 'Blues', 'Grens', 'Purples', "Oranges", 'Greys']

over9 = {'Paired', 'Paired_r', 'Set3', 'Set3_r'}
over8 = over9 | {'Set1', "Pastel1"}

def extract_label_array(table, tax_col, tax_delim='|'):
    """
    Converts delimited taxonomy strings into a working table
    
    Parameters
    ----------
    table : DataFrame
        A DataFrame with observation on the rows (biom-style table) with 
        `tax_col` as one of its columns.
    tax_col : str
        The column in `table` containig the taxobnomy information
    tax_delim: str, optional
        The delimiter between taxonomic groups
        
    Returns
    -------
    DataFrame
        The taxonomic strings parsed into n levels
    """
    def f_(x):
        return pd.Series([y.strip() for y in x.split(tax_delim)])
    
    return table[tax_col].apply(f_)


def level_taxonomy(table, taxa, samples, level, consider_nan=True):
    """
    Gets the taxonomy collapsed to the desired level
    
    Parameters
    ----------
    table : DataFrame
        A table with observation on the rows (biom-style table) with `samples`
        in its columns.
    taxa: DataFrame
        The taxonomic strings parsed into n levels
    level: int
        The level to which the taxonomy should be summarized 
    samples : list
        The columns from `table` to be included in the analysis
    consider_nan: bool, optional
        Whether the table contains multiple concatenated, in which cases 
        considering `nan` will filter the samples to retain only the levels 
        of interest. This is recommended for kracken/bracken tables, but not 
        applicable for some 16s sequences
    """
    level = level.max()
    if consider_nan:
        leveler = (taxa[level].notna() & taxa[(level + 1)].isna())
    else:
        leveler = (taxa[level].notna())
            
    cols = list(np.arange(level + 1))
    
    # Combines the filtered tables
    level_ = pd.concat(
        axis=1, 
        objs=[taxa.loc[leveler, cols],
              (table.loc[leveler, samples] / 
               table.loc[leveler, samples].sum(axis=0))],
        # sort=False
    )
    if taxa.loc[leveler, cols].duplicated().any():
        return level_.groupby(cols).sum()
    else:
        return level_.set_index(cols)


def profile_one_level(collapsed, level, threshhold=0.01, count=8):
    """
    Gets upper and lower tables for a single taxonmic level
    
    Parameters
    ----------
    Collapsed: DataFrame
        The counts data with the index as a multi-level index of 
        levels of interest and the columns as samples
    threshhold: float, optional
        The minimum relative abundance for an organism to be shown
    count : int, optional
        The maximum number of levels to show for a single group
        
    Returns
    -------
    DataFame
        A table of the top taxa for the data of interest
    """
    collapsed['mean'] = collapsed.mean(axis=1)
    collapsed.sort_values(['mean'], ascending=False, inplace=True)
    collapsed['count'] = 1
    collapsed['count'] = collapsed['count'].cumsum()

    thresh_ = (collapsed['mean'] > threshhold) & (collapsed['count'] <= count)
    top_taxa = collapsed.loc[thresh_].copy()
    top_taxa.drop(columns=['mean', 'count'], inplace=True)
    
    top_taxa.reset_index(inplace=True)
    top_taxa.drop(columns=np.arange(level), inplace=True)
    top_taxa.set_index(level, inplace=True)
    
    top_taxa = top_taxa[samples].sort_values(by=[top_taxa.index[0]], 
                                             ascending=False,
                                             axis='columns')

    upper_ = top_taxa.cumsum()
    lower_ = top_taxa.cumsum() - top_taxa

    return upper_, lower_


def profile_joint_levels(collapsed, lo_, hi_, lo_thresh=0.01, 
    hi_thresh=0.01, lo_count=4, hi_count=5):
    """
    Generates a table of taxonomy using two levels to define grouping
    
    Parameters
    ----------
    collapsed: DataFrame
        The counts data with the index as a multi-level index of 
        levels of interest and the columns as samples
    lo_, hi_: int
        The numeric identifier for lower (`lo_`) and higher (`hi_`)
        resolution where "low" is defined as having fewer groups.
        (i.e. for taxonomy Phylum is low, family is high)
    lo_thresh, hi_thresh: int, optional
        The minimum relative abundance for an organism to be shown
        at a given level
    lo_count, hi_count : int, optional
        The maximum number of levels to show for a single group. 
        This is to appease the limitations of our eyes and colormaps.
        
    Returns
    -------
    DataFame
        A table of the top taxa for the data of interest
    """
    collapsed['mean_hi'] = collapsed.mean(axis=1)
    collapsed.reset_index(inplace=True)
    mean_lo_rep = collapsed.groupby(lo_)['mean_hi'].sum().to_dict()
    collapsed['mean_lo'] = collapsed[lo_].replace(mean_lo_rep)
    collapsed.sort_values(['mean_lo', 'mean_hi'], ascending=False, inplace=True)

    collapsed['count_lo'] = ~collapsed[lo_].duplicated(keep='first') * 1.
    collapsed['count_lo'] = collapsed['count_lo'].cumsum()
    collapsed['count_hi'] = 1
    collapsed['count_hi'] = collapsed.groupby(lo_)['count_hi'].cumsum()

    collapsed['thresh_lo'] = ((collapsed['mean_lo'] > lo_thresh) & 
                              (collapsed['count_lo'] <= lo_count))
    collapsed['thresh_hi'] = ((collapsed['mean_hi'] > hi_thresh) & 
                              (collapsed['count_hi'] <= hi_count))

    top_lo = collapsed.loc[collapsed['thresh_lo']].copy()
    top_lo['other'] = ~(top_lo['thresh_lo'] & ~top_lo['thresh_hi']) * 1
    top_lo['new_name'] = top_lo[lo_].apply(lambda x: 'other %s' % x )
    top_lo.loc[top_lo['thresh_hi'], 'new_name'] = \
        top_lo.loc[top_lo['thresh_hi'], hi_]

    new_taxa = top_lo.groupby([lo_, 'new_name']).sum()
    new_taxa.sort_values(by=['mean_lo', 'other', 'mean_hi'], 
                         ascending=False, inplace=True)
    
    upper_ = new_taxa.cumsum()
    lower_ = new_taxa.cumsum() - new_taxa

    return upper_, lower_


def define_single_cmap(cmap, top_taxa):
    """
    Gets the colormap a single level table

    """
    # Gets the colormap object
    map_ = cm.cmap_d[cmap]
    # Gets the taxonomic object
    return {tax: map_(i) for i, tax in enumerate(top_taxa.index)}


def define_join_cmap(table):
    """
    Defines a joint colormap for a taxonomic table.
    """
    grouping = table['dummy'].reset_index()

    rough_order = grouping['rough'].unique()

    rough_map = {group: cm.cmap_d[cmap] 
                 for (group,cmap) in zip(*(rough_order, colors_order))} 
    pooled_map = dict([])
    for rough_, fine_ in grouping.groupby('rough')['fine']:
        cmap_ = rough_map[rough_]
        colors = {c: cmap_(200-(i + 1) * 20) for i, c in enumerate(fine_)}
        pooled_map.update(colors)

    return pooled_map


def plot_area(upper_, lower_, colors, sample_interval=5):
    """
    An in-elegant function to make an area plot

    Yes, you'll get far more control if you do it yourself outside this
    function but it will at least give you a first pass of a stacked area
    plot
    
    Parameters
    ---------
    upper_, lower_ : DataFrame
        The upper (`top_`) and lower (`low_`) limits for the 
        area plot. This should already be sorted in the
        desired order.
    colors: dict
        A dictionary mapping the taxonomic label to the
        appropriate matplotlib readable colors. For convenience,
        `define_single_colormap` and `define_joint_colormap`
        are good functions to use to generate this
    sample_interval : int, optional
        The interval for ticks for counting samples.
        
    Returns
    -------
    Figure
        A 8" x 4" matplotlib figure with the area plot and legend.
    """
    
    # Gets the figure
    fig_, ax1 = plt.subplots(1,1)
    fig_.set_size_inches((8, 4))
    ax1.set_position((0.15, 0.125, 0.4, 0.75))
    
    # Pl;ots the area plot
    x = np.arange(0, len(top_.columns))
    for taxa, hi_ in upper_.iterrows():
        lo_ = lower_.loc[taxa]
        cl_ = colors[taxa]
    
        ax1.fill_between(x=x, y1=1-lo_, y2=1-hi_, 
                         color=cl_, label=taxa)
    # Adds the legend
    leg_ = ax1.legend()
    leg_.set_bbox_to_anchor((2.05, 1))

    # Sets up the y-axis so the order matches the colormap
    # (accomplished by flipping the axis?)
    ax1.set_ylim((1, 0))
    ax1.set_yticks(np.arange(0, 1.1, 0.25))
    ax1.set_yticklabels(np.arange(1, -0.1, -0.25), size=11)
    ax1.set_ylabel('Relative Abundance', size=13)

    # Sets up x-axis without numeric labels
    ax1.set_xticklabels([])
    ax1.set_xticks(np.arange(0, len(samples), sample_interval))
    ax1.set_xlim((0, len(samples) - 1))
    ax1.set_xlabel('Samples', size=13)

    return fig_


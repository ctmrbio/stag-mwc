#!/usr/bin/env python3

from sys import argv, exit
import argparse
import warnings

from matplotlib import rcParams
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

"""
Generates a pretty areaplot from a collapsed feature table.
"""

__author__ = 'JW Debelius'
__date__ = '2020-02'
__version__ = "0.2"

# Sets up the matplotlib parameters so that we can save to be edited in 
# illustator if a direct conversion is required. Because it just makes life
# better
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Sets up the order of colors to be used in joint plots.
colors_order = ['Reds', 'Blues', 'Greens', 'Purples', "Oranges", 'Greys']

over9 = {'Paired', 'Paired_r', 'Set3', 'Set3_r'}
over8 = over9 | {'Set1', "Pastel1"}

mode_dict = {
    'kaiju': {
        'tax_delim': ';',
        'multi_level': True,
        'tax_col': 'taxon_name',
        'table_drop': [],
        'skip_rows': 0,
    },
    'kraken2': {
        'tax_delim': '|',
        'multi_level': True,
        'tax_col': 'taxon_name',
        'table_drop': [],
        'skip_rows': 0,
    },
    'metaphlan': {
        'tax_delim': '|',
        'multi_level': True,
        'tax_col': 'clade_name',
        'table_drop': ['NCBI_tax_id'],
        'skip_rows': 1,
    },
    'metaphlan4': {
        'tax_delim': '|',
        'multi_level': True,
        'tax_col': 'clade_name',
        'table_drop': [],
        'skip_rows': 1,
    },
    'marker': {
        'tax_delim': ';',
        'multi_level': False,
        'tax_col': 'taxonomy',
        'table_drop': ['sequence'],
        'skip_rows': 0,
    },
}


def extract_label_array(table, tax_col, tax_delim='|'):
    """
    Converts delimited taxonomy strings into a working table
    
    Parameters
    ----------
    table : DataFrame
        A DataFrame with observation on the rows (biom-style table) with 
        `tax_col` as one of its columns.
    tax_col : str
        The column in `table` containing the taxonomy information
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
    level: list
        The level to which the taxonomy should be summarized 
    samples : list
        The columns from `table` to be included in the analysis
    consider_nan: bool, optional
        Whether the table contains multiple concatenated, in which cases 
        considering `nan` will filter the samples to retain only the levels 
        of interest. This is recommended for kraken/bracken tables, but not 
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
              (table.loc[leveler, samples] / table.loc[leveler, samples].sum(axis=0))],
    #     # sort=False
    )
    level_.reset_index()
    if taxa.loc[leveler, cols].duplicated().any():
        return level_.groupby(cols).sum()
    else:
        return level_.set_index(cols)


def profile_one_level(collapsed, level, threshold=0.01, count=8):
    """
    Gets upper and lower tables for a single taxonomic level
    
    Parameters
    ----------
    Collapsed: DataFrame
        The counts data with the index as a multi-level index of 
        levels of interest and the columns as samples
    threshold: float, optional
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

    thresh_ = (collapsed['mean'] > threshold) & (collapsed['count'] <= count)
    top_taxa = collapsed.loc[thresh_].copy()
    top_taxa.drop(columns=['mean', 'count'], inplace=True)
    for l_ in np.arange(level):
        top_taxa.index = top_taxa.index.droplevel(l_)

    first_ = top_taxa.index[0]
    
    top_taxa.sort_values(
        [first_],
        ascending=False,
        axis='columns',
        inplace=True,
    )

    upper_ = top_taxa.cumsum()
    lower_ = top_taxa.cumsum() - top_taxa

    return upper_, lower_


def profile_joint_levels(collapsed, lo_, hi_, samples, lo_thresh=0.01, 
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

    drop_levels = np.arange(hi_)[np.arange(hi_) != lo_]
    top_lo.drop(columns=drop_levels, inplace=True)
    new_taxa = top_lo.groupby([lo_, 'new_name']).sum(dropna=True)
    new_taxa.reset_index(inplace=True)
    new_taxa['mean_lo'] = new_taxa[lo_].replace(mean_lo_rep)
    new_taxa.set_index([lo_, 'new_name'], inplace=True)
    new_taxa.sort_values(['mean_lo', 'count_hi', 'mean_hi'], 
                         ascending=[False, True, False], 
                         inplace=True)

    upper_ = new_taxa.cumsum()[samples]
    upper_.sort_values([upper_.index[0], upper_.index[1]],
                       axis='columns', 
                       inplace=True, ascending=False)
    lower_ = upper_ - new_taxa[upper_.columns]
    
    upper_.index.set_names(['rough', 'fine'], inplace=True)
    lower_.index.set_names(['rough', 'fine'], inplace=True)

    return upper_, lower_


def define_single_cmap(cmap, top_taxa):
    """
    Gets the colormap a single level table

    """
    # Gets the colormap object
    map_ = mpl.colormaps[cmap]
    # Gets the taxonomic object
    return {tax: map_(i) for i, tax in enumerate(top_taxa.index)}


def define_join_cmap(table):
    """
    Defines a joint colormap for a taxonomic table.
    """
    table['dummy'] = 1
    grouping = table['dummy'].reset_index()

    rough_order = grouping['rough'].unique()

    rough_map = {group: mpl.colormaps[cmap] 
                 for (group,cmap) in zip(*(rough_order, colors_order))} 
    pooled_map = dict([])
    for rough_, fine_ in grouping.groupby('rough')['fine']:
        cmap_ = rough_map[rough_]
        colors = {c: cmap_(200-(i + 1) * 20) for i, c in enumerate(fine_)}
        pooled_map.update(colors)

    table.drop(columns=['dummy'], inplace=True)

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
    
    # Plots the area plot
    x = np.arange(0, len(upper_.columns))
    for taxa, hi_ in upper_.iloc[::-1].iterrows():
        lo_ = lower_.loc[taxa]
        cl_ = colors[taxa]

    
        ax1.fill_between(x=x, y1=1-lo_.values, y2=1-hi_.values, 
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
    ax1.set_xticks(np.arange(0, x.max(), sample_interval))
    ax1.set_xlim((0, x.max() - 0.99))  # Subtract less than 1 to avoid singularity if xmin=xmax=0
    ax1.set_xlabel('Samples', size=13)

    return fig_


def single_area_plot(table, level=3, samples=None, 
    tax_col='taxon_name', cmap='Set3',
    tax_delim='|', multilevel_table=True, abund_thresh=0.1, 
    group_thresh=8):
    """
    Generates an area plot for the table at the specified level of resolution

    Parameters
    ----------
    table : DataFrame
        A pandas dataframe of the original table of data (either containing 
        counts or relative abundance)
    level : int
        The hierarchical level within the table to display as an integer
    cmap : str
        The qualitative colormap to use to generate your plot. Refer to 
        colorbrewer for options. If a selected colormap exceeds the number
        of groups (`--group-thresh`) possible, it will default to Set3.
    samples : list, optional
        The columns from `table` to be included in the analysis. If `samples`
        is None, then all columns in `table` except `tax_col` will be used.
    tax_col : str, optional
        The column in `table` which contains the taxonomic information.
    tax_delim: str, optional
        The delimiter between taxonomic levels, for example "|" or ";".
    multilevel_table: bool, optional
        Whether the table contains multiple concatenated, in which cases 
        considering `nan` will filter the samples to retain only the levels 
        of interest. This is recommended for kraken/bracken tables, but not 
        applicable for some 16s sequences
    abund_thresh: float [0, 1]
        The mean abundance threshold for a sample to be plotted. This is 
        in conjunction with the group threshold (`--group-thresh`) will be 
        used to determine the groups that are shown.
    group_thresh: int, [1, 12]
        The maximum number of groups (colors) to show in the area plot. This 
        is handled in conjunction with the `--abund-thresh` in that 
        to be displayed, a group must have both a mean relative abundance 
        exceeding the `abund-thresh` and must be in the top `group-thresh` 
        groups.

    Returns
    -------
    Figure
        A 8" x 4" matplotlib figure with the area plot and legend.

    Also See
    --------
    make_joint_area_plot
    
    """

    if group_thresh > 12:
        raise ValueError("You may display at most 12 colors on this plot. "
                         "Please re-consider your plotting choices.")
    elif (group_thresh > 9) & ~(cmap in over9):
        warnings.warn('There are too many colors for your colormap. '
                      'Changing to Set3.')
        cmap = 'Set3'
    elif (group_thresh > 8) & ~(cmap in over8):
        warnings.warn('There are too many colors for your colormap. '
                      'Changing to Set3.')
        cmap = 'Set3'

    # Parses the taxonomy and collapses the table
    taxa = extract_label_array(table, tax_col, tax_delim)

    if samples is None:
        samples = list(table.columns.values)
        samples.remove(tax_col)

    # Gets the appropriate taxonomic level information to go forward
    collapsed = level_taxonomy(table, taxa, samples, np.array([level]), 
                              consider_nan=multilevel_table)

    # Gets the top taxonomic levels
    upper_, lower_, = profile_one_level(collapsed, np.array([level]), 
                                        threshold=abund_thresh, 
                                        count=group_thresh)

    # Gets the colormap 
    cmap = define_single_cmap(cmap, upper_)

    # Plots the data
    fig_ = plot_area(upper_, lower_, cmap)

    return fig_


def joint_area_plot(table, rough_level=2, fine_level=5, samples=None, 
    tax_col='taxon_name', tax_delim='|', 
    multilevel_table=True, abund_thresh_rough=0.1, 
    abund_thresh_fine=0.05, group_thresh_fine=5, 
    group_thresh_rough=5):
    """
    Generates an area plot with nested grouping where the the higher level
    (`rough_level`) in the table (lower resolution/fewer groups) is used to 
    provide the general grouping structure and then within each `rough_level`,
    a number of `fine_level` groups are displayed. 

    Parameters
    ----------
    table : DataFrame
        A dataframe of hte original data, either as counts or relative 
        abundance with the taxonomic information in `tax_col`. The data
        can have separate count values at multiple levels (i.e. combine)
        collapsed phylum, class, etc levels.
    rough_level, fine_level: int
        The taxonomic levels to be displayed. The `fine_level` will be grouped
        by `rough_level` to display the data grouped by `rough_level`. The
        `rough_level` should smaller than the `fine_level`. 
    samples : list, optional
        The columns from `table` to be included in the analysis. If `samples`
        is None, then all columns in `table` except `tax_col` will be used.
    tax_col : str, optional
        The column in `table` which contains the taxonomic information.
    tax_delim: str, optional
        The delimiter between taxonomic levels, for example "|" or ";".
    multilevel_table: bool, optional
        Whether the table contains multiple concatenated, in which cases 
        considering `nan` will filter the samples to retain only the levels 
        of interest. This is recommended for kraken/bracken tables, but not 
        applicable for some 16s sequences
    abund_thresh_rough, abund_thresh_fine : float [0, 1]
        The mean abundance threshold for a taxonomic group to be plotted for
        the higher level grouping (`abund_thresh_rough`) and sub grouping
        level. This will be used in conjunction with the `group_thresh_rough`
        and `group_thresh_fine` to determine the number of groups to be
        included.
    group_thresh_fine, group_thresh_rough: int, [1, 6]
        The maximum number of taxonmic groups to display for the respective 
        level. If `group_thresh_rough` > 6, then it will be replaced with 
        6 because this is the maximum number of avaliable color groups.

    Returns
    -------
    Figure
        A 8" x 4" matplotlib figure with the area plot and legend.

    Also See
    --------
    single_area_plot

    """

    # Parses the taxonomy and collapses the table
    taxa = extract_label_array(table, tax_col, tax_delim)
    if samples is None:
        samples = list(table.drop(columns=[tax_col]).columns)

    # Gets the appropriate taxonomic level information to go forward
    collapsed = level_taxonomy(table, taxa, samples, 
                               level=np.array([fine_level]), 
                               consider_nan=multilevel_table)
    samples = collapsed.columns

    # Gets the top taxonomic levels
    upper_, lower_, = profile_joint_levels(collapsed, rough_level, fine_level, 
                                           samples=samples,
                                           lo_thresh=abund_thresh_rough, 
                                           lo_count=min(5, group_thresh_rough),
                                           hi_thresh=abund_thresh_fine,
                                           hi_count=group_thresh_fine,
                                           )
    # Gets the colormap 
    cmap = define_join_cmap(upper_)
    upper_.index = upper_.index.droplevel('rough')
    lower_.index = lower_.index.droplevel('rough')

    # Plots the data
    fig_ = plot_area(upper_.astype(float), lower_.astype(float), cmap)

    return fig_


# Sets up the main arguments for argparse.
def create_argparse():
    parser_one = argparse.ArgumentParser(
        description=('A set of functions to generate diagnostic stacked area '
                     'plots from metagenomic outputs.'),
        prog=('area_plotter'),
        )
    parser_one.add_argument(
        '-t', '--table', 
        help=('The abundance table as a tsv classic biom (features as rows, '
              'samples as columns) containing absloute or relative abundance '
              'for the samples.'),
        required=True,
        )
    parser_one.add_argument(
        '-o', '--output',
        help=('The location for the final figure'),
        required=True,
        )
    parser_one.add_argument(
        '-s', '--samples', 
        help=('A text file with the list of samples to be included (one '
            'per line). If no list is provided, then data from all columns '
            'in the table (except the one specifying taxonomy) will be used.'),
        )
    parser_one.add_argument(
        '--mode', 
        choices=mode_dict.keys(),
        help=('The software generating the table to make parsing easier. '
              'Options are kraken, metaphlan, marker (i.e. CTMR amplicon).'),
        )
    parser_one.add_argument(
        '-l', '--level',
        help=('The taxonomic level (as an integer) to plot the data.'),
        default=3,
        type=int,
        )
    parser_one.add_argument(
        '--abund-thresh',
        help=("the minimum abundance required to display a group."),
        default=0.01,
        type=float,
        )
    parser_one.add_argument(
        '--group-thresh',
        help=("The maximum number of groups to be displayed in the graph."),
        default=8,
        type=int,
        )
    parser_one.add_argument(
        '-c', '--colormap',
        help=("The qualitative colormap to use to generate your plot. Refer"
             ' to colorbrewer for options. If a selected colormap exceeds '
             'the number of groups (`--group-thresh`) possible, it will '
             'default to Set3.'),
        default='Set3',
        )
    parser_one.add_argument(
        '--sub-level',
        help=('The second level to use if doing a joint plot'),
        type=int,
        )
    parser_one.add_argument(
        '--sub-abund-thresh',
        help=("the minimum abundance required to display a sub group"),
        default=0.05,
        type=float,
        )
    parser_one.add_argument(
        '--sub-group-thresh',
        help=("the maximum number of sub groups allowed in a joint level plot."),
        default=5,
        type=float,
        )
    parser_one.add_argument(
        '--tax-delim',
        help=("String delimiting taxonomic levels."),
        type=str,
        )
    parser_one.add_argument(
        '--multi-level',
        help=("Whether the table contains multiple concatenated, in which "
             "case considering `nan` will filter the samples to retain only"
             "the levels of interest. This is recommended for most "
             "metagenomic tables, but not applicable for some 16s sequences."),
        )
    parser_one.add_argument(
        "--tax-col",
        help=("The column in `table` containig the taxobnomy information"),
        )
    parser_one.add_argument(
        '--table-drop',
        help=('A comma-seperated list describing the columns to drop'),
        )
    parser_one.add_argument(
        '--skip-rows',
        help=('The number of rows to skip when reading in the feature table.')
        )

    return parser_one


if __name__ == '__main__':
    parser_one = create_argparse()

    if len(argv) < 2:
        parser_one.print_help()
        exit()

    args = parser_one.parse_args()

    if args.table_drop is not None:
        args.table_drop = [s for s in args.table_drop.split(',')]
    else:
        args.table_drop = []

    mode_defaults = mode_dict.get(args.mode, mode_dict['kraken2'])
    mode_defaults.update({k: v for k, v in args.__dict__.items() 
                         if (k in mode_defaults) and (v)})

    table = pd.read_csv(args.table, sep='\t', 
                        skiprows=mode_defaults['skip_rows'])
    if args.samples is not None:
        with open(args.samples, 'r') as f_:
            samples = f_.read().split('\n')
    else:
        samples = None

    if args.sub_level is not None:
        fig_ = joint_area_plot(
            table.drop(columns=mode_defaults['table_drop']),
            rough_level=args.level,
            fine_level=args.sub_level,
            samples=args.samples,
            tax_delim=mode_defaults['tax_delim'],
            tax_col=mode_defaults['tax_col'],
            multilevel_table=mode_defaults['multi_level'],
            abund_thresh_rough=args.abund_thresh,
            group_thresh_rough=args.group_thresh,
            abund_thresh_fine=args.sub_abund_thresh,
            group_thresh_fine=args.sub_group_thresh,
        )
    else:
        fig_ = single_area_plot(
            table.drop(columns=mode_defaults['table_drop']),
            level=args.level,
            cmap=args.colormap,
            samples=samples,
            tax_delim=mode_defaults['tax_delim'],
            tax_col=mode_defaults['tax_col'],
            multilevel_table=mode_defaults['multi_level'],
            abund_thresh=args.abund_thresh,
            group_thresh=args.group_thresh,
        )

    fig_.savefig(args.output, dpi=300)

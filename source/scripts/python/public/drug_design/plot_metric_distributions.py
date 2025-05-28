#!/usr/bin/env python2.7
##
## (c) Copyright Rosetta Commons Member Institutions.
## (c) This file is part of the Rosetta software suite and is made available under license.
## (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
## (c) For more information, see http://www.rosettacommons.org. Questions about this can be
## (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
##
## @author Rocco Moretti (rmorettiase@gmail.com)

'''
Plots the distribution of metrics of experimental values compared to other distributions.

Default color order is: blue, red, orange, yellow, green. (Unneeded colors won't be used.)

Output can be in either pdf or png format (chosen by output extension). If PDF, all plots are in the same file.
If PNG, a series of files will be made, each with filenames based on the given name expanded with the metric name.

While running, it will print which values it's plotting to the standard output.

Requires Pandas
'''

import argparse
import pandas
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys, os

pandas.set_option('mode.use_inf_as_null', True)

lighten_color_dict = {
    "blue":"#3030FF",
    "red":'#FF3030',
    "orange":'#FFB220',
    "yellow":"#FFFF30",
    "green":"#159015",
    "b":"#3030FF",
    "r":'#FF3030',
    "o":'#FFB220',
    "y":"#FFFF30",
    "g":"#159015",
}

color_trans_dict = {
    "o": "orange",
}

def lighten_color(color):
    return lighten_color_dict.get(color,color)

def get_colors(input_colors):
    if input_colors is None:
        return ["blue","red","orange","yellow","green"]
    colors = input_colors.split(',')
    return [ color_trans_dict.get(c,c) for c in colors ]

def change_column_names( frame ):
   replacements = [ ( 'MoleculeSum(IsNotH)', "NHeavy" ),
     ( 'Subtract(lhs=MoleculeSum(Abs(Atom_Stereocenters)),rhs=MoleculeSum(Equal(Atom_Stereocenters,2)))', "NChiral" ),
     ( 'TopologicalPolarSurfaceArea', 'TPSA' ),
   ]
   cn = list(frame.columns)
   for val, rep in replacements:
      if val in cn:
       cn[ cn.index( val ) ] = rep
   frame.columns = cn

def do_plot(value, all):
    all_vec = all[value]
    # Don't plot non-numeric values
    if type(all_vec.max()) == str or type(all_vec.min()) == str: return False
    # Don't plot values with really large values: something is probably wrong with the calculation
    if all_vec.max() > 1e20: return False

    return True

def do_decorrelate(value, all, decorrelation=None):
    if decorrelation is None: return False
    if value == decorrelation: return False
    if decorrelation not in all:
        print "ERROR: Can't find metric '%s' for decorrelation." % str( decorrelation )
        return False

    cor = all[value].corr(all[decorrelation])
    if abs(cor) <= 0.5: return False
    cor2 = (all[value]/all[decorrelation]).corr(all[decorrelation])
    if abs(cor2) >= abs(cor): return False # Don't do 'decorrelation' if the decorrelated values are actually better correlated.

    return True

def plot_metric_figure(value, all, tables, skip_outlier_detect, quantile, decorrelation=None, colors=None, labels=None):
    if colors is None:
        colors = get_colors(None)

    if not do_plot(value, all): return True
    decorrelate = do_decorrelate(value, all, decorrelation)
    if decorrelation is not None and not decorrelate: return True # If we're trying to decorrelate, but we shouldn't/can't

    all_vec = all[value]

    if decorrelate:
        all_vec = all_vec / all[decorrelation]

    all_vec = all_vec[ all_vec.notnull() ]

    # Don't plot values with no variance across any entry (automatic range finding doesn't work)
    if( all_vec.max() - all_vec.min() == 0 ):
        print "Skipping", value
        return True
    # Remove large outliers
    sd = all_vec.std()
    m = all_vec.mean()
    if not skip_outlier_detect:
        all_vec = all_vec[ all_vec < m + 3 * sd ]
        all_vec = all_vec[ all_vec > m - 3 * sd ]
    if 0 < quantile < 1.0:
        trim = (1.0-quantile)/2.0
        lower = all_vec.quantile(trim)
        upper = all_vec.quantile(1.0-trim)
        all_vec = all_vec[ all_vec >= lower ]
        all_vec = all_vec[ all_vec <= upper ]
    elif quantile != 1.0:
        raise ValueError("Quantile value " + str( quantile ) + " is not a valid value to use.")

    unit = (all_vec.max() - all_vec.min())/100.0
    if( unit == 0 ):
        print "Skipping", value
        return True

    label = value
    if decorrelate:
        label = value + " / " + decorrelation

    print label

    value_series = []
    for source in tables:
        vec = source[value]
        if decorrelate:
            vec = vec / source[decorrelation]
        vec = vec[ vec.notnull() ]
        vec = vec[ vec >= all_vec.min() ]
        vec = vec[ vec <= all_vec.max() ]
        value_series.append(vec)

    if all_vec.dtype == int and  all_vec.min() >= 0 and all_vec.max() <= 10:
        # Small number of unique values? Use (normed) histogram plots:
        colors = [ lighten_color(c) for c in colors ][:len(tables)]
        plt.hist(value_series, bins=[ x-0.5 for x in range(all_vec.max()+2) ], color=colors, normed=True, align="mid", rwidth=0.67 )
        plt.xticks( range(all_vec.max()+1) )
    else:
        # Large number of unique values? Use dist plots.

        # Find the global bandwith for consistentcy
        all_kde = scipy.stats.gaussian_kde(all_vec)
        bw = all_kde.factor*all_vec.std() # Dumb SciPy does factor instead of actual bandwidth ...
        xlimits = [all_vec.min() - 3*bw, all_vec.max() + 3*bw ]

        ind = np.linspace(xlimits[0], xlimits[1], 1000)

        for j, vec in enumerate(value_series):
            if( vec.std() > unit ):
                kde = scipy.stats.gaussian_kde(vec, bw_method=bw/vec.std())
                plt.plot(ind, kde.evaluate(ind), color=colors[j] )
            else:
                plt.axvline( vec.mean(), color=colors[j])
        plt.gca().set_xlim( xlimits )
    plt.title(label)
    plt.xlabel(label)
    if labels is not None:
        plt.legend( labels, loc=0 )
    return False

def main(args):
    outname = args.output_file
    ext = os.path.splitext(outname)[1]
    pdfname = ""
    pngstem = ""
    if ext.lower() == ".pdf":
        pdfname = outname
    elif ext.lower() == ".png":
        pngstem = os.path.splitext(outname)[0]
    else:
        raise ValueError("Can't understand output format for '" + outname + "' - is not pdf/png.")

    tables = [ pandas.read_table(f, sep="[\s]+", engine='python' ) for f in args.tables ]

    colors = get_colors(args.colors)
    if len(tables) > len(colors):
       raise ValueError("Color vector too short for the number of tables. (" + str(len(tables)) + " tables but only " + str(len(colors)) + " colors.)" )

    labels = args.labels
    if labels is not None:
        labels = args.labels.split(';')
        if len(labels) != len(tables):
            raise ValueError("Labels length must match the number of tables used. (" + str(len(tables)) + " tables but " + str(len(labels)) + " labels.)" )

    for table in tables:
        change_column_names( table )

    all = pandas.concat( tables )

    metrics = (args.metrics or [])
    if args.list is not None:
        with open(args.list) as f:
            for line in f:
                line = line.strip()
                if len(line) > 0 and line[0] != '#':
                    metrics.append( line )
    for m in metrics:
        if m not in all.columns:
            raise ValueError("Metric '"+m+"' is not present in the input files!")

    if len(metrics) == 0:
        metrics = all.columns

    if pdfname:
        pp = PdfPages( pdfname )

    i = 1
    for metric in metrics:
        if metric == "Index":
            continue

        if not do_plot(metric, all):
            continue

        if args.decorrelate is None:
            decorrs = [None]
        elif args.both:
            decorrs = [None, args.decorrelate]
        else:
            decorrs = [args.decorrelate]

        for decorr in decorrs:
            if decorr is not None and not do_decorrelate(metric, all, decorr): continue

            plt.figure(i)

            if plot_metric_figure(metric, all, tables, args.no_outlier, args.q, decorr, colors, labels):
                continue

            if pdfname:
                plt.savefig(pp, format='pdf')
            elif pngstem:
                plt.savefig(pngstem + '_' + metric + '.png' )
            plt.close(i)
            i += 1

    if pdfname:
        pp.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('output_file', help="The output file to use. Use extension pdf/png to change format.")
    parser.add_argument('tables', nargs='*', help="The metric tables to plot.")
    parser.add_argument('--colors', default=None, help="A comma-separated list of colors to use.")
    parser.add_argument('-l', '--labels', default=None, help="A semicolon-separated list of labels to use for the plots.")
    parser.add_argument('--no_outlier', '-o', action='store_true', help="Skip outlier removal when determining plotting range. (Outliers are > 3sd outside)")
    parser.add_argument('-q', type=float, default=1.0, help="Scale the output window for the inner q quantile of values.")
    parser.add_argument('--decorrelate', help="Divide through by this value (e.g. MolWt) before plotting if the metric is highly correlated.")
    parser.add_argument('-b','--both', action='store_true', help="When decorrelating, plot both the decorrelated and the non-decorrelated values.")
    parser.add_argument('--metrics', nargs='*', help="Which metrics to plot (defaults to all). If using, place this last.")
    parser.add_argument('--list', help="The name of a file which contains a list of metrics to plot")

    args = parser.parse_args()
    main(args)

#!/usr/bin/env python

# This code was well-written, but it's starting to break down.  It does so much 
# more than it was designed to.  

from __future__ import division

import sys, os, math, pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(
            description="""\
                    Display a histogram showing which torsion angles were 
                    sampled during a particular kinematic closure trajectory.  
                    Arguments are provided to control how much data is read off 
                    of disk (the slowest step), how the data is analyzed, and 
                    how the analysis is displayed. """,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('jobs', nargs='+',
            help="Directory containing the data files.")

    parser.add_argument('--include', '-i', nargs='+', default=[],
            help="Show only the given torsions " + '(e.g phi0, psi3, omega).')

    parser.add_argument('--exclude', '-x', nargs='+', default=['omega3'],
            help="Hide the given torsions " + '(e.g phi0, psi3, omega).')

    parser.add_argument('--force-refresh', '-f', action='store_true',
            help="Force the numpy torsion cache to be refreshed.")
    
    parser.add_argument('--bins', '-b', type=int, default=360,
            help="Number of bins to histogram the data on.")
    
    parser.add_argument('--normalize', '-n', action='store_true',
            help="Normalize the histograms before plotting.")
    
    parser.add_argument('--limit', '-l', type=int,
            help="Maximum number of iterations to read off of disk.")
    
    parser.add_argument('--title', '-t',
            help="Title to put on the resulting histogram.")

    parser.add_argument('--ylim', '-y', nargs=2,
            help="Minimum and maximum values of the vertical axis.")
    
    parser.add_argument('--ks-test', '-ks', action='store_true',
            help="Show the KS statistic for each histogram.")

    parser.add_argument('--rama', '-r', action='store_true',
            help="Compare the trajectories to the Ramachandran distribution.")

    parser.add_argument('--style', '-s', nargs='+', default=['rainbow'],
            help="Specify the style that should be used for each histogram.")
    
    parser.add_argument('--output', '-o',
            help="File path where the histogram should be saved.")
    
    parser.add_argument('--dpi', type=int, default=300,
            help="Resolution to use for raster output formats.")
    
    parser.add_argument('--no-gui', '-q', action='store_true',
            help="Indicate that a GUI histogram should not be displayed.")
    
    arguments = parser.parse_args()

    if arguments.rama:
        arguments.normalize = True

    return arguments

def parse_iterations(arguments):
    """ Return the number of iterations that need to be read off disk. """
    job = os.path.join(arguments.job, 'job.sh')

    if arguments.limit:
        return arguments.limit
    else:
        with open(job) as file:
            for line in file:
                fields = line.strip().split()
                if fields and fields[0] == '-kale:mc:iterations':
                    return fields[1]

def key_from_torsion(torsion, index):
    """ Return a simple string representing the given torsion angle.  A handful 
    of command line options expect arguments in this format. """
    return torsion + str(index)

def torsion_from_key(key):
    """ Return a (torsion, index) tuple represented by this key. """
    import re
    return re.match('(\w+)(\d+)', key).groups()


def load_histograms(arguments):
    """ Return two dictionaries holding the bins and the counts, respectively, 
    for each torsion.  Much effort has been made to optimize this function by 
    caching intermediate results. """

    directory = os.path.join(arguments.job, 'histogram')
    try: os.makedirs(directory)
    except OSError:
        if os.path.isdir(directory): pass
        else: raise
    
    # The user can disable all caching by providing the '--force-refresh' flag.  
    # If this flag has been given, just regenerate the histograms and don't 
    # bother with all the logic below.

    if arguments.force_refresh:
        torsions = load_torsions_from_dat(arguments)
        return generate_histograms(arguments, torsions)

    # First try to load the histograms themselves directly from disk.  This is 
    # extremely fast, because the histograms themselves are quite tiny.

    try:
        return load_histograms_from_npz(arguments)
    except IOError as error:
        pass

    # If the histograms themselves are not present, try loading the compressed 
    # data array before reading the data from a text file.  If the compressed 
    # data is present, it will load much faster than the text file.

    try:
        torsions = load_torsions_from_npz(arguments)
    except IOError as error:
        torsions = load_torsions_from_dat(arguments)

    return generate_histograms(arguments, torsions)

def load_histograms_from_npz(arguments):
    bin_path = os.path.join(arguments.job, 'histogram', 'bins.npz')
    count_path = os.path.join(arguments.job, 'histogram', 'counts.npz')
    ks_path = os.path.join(arguments.job, 'histogram', 'ks_test.pkl')

    with open(ks_path) as file:
        ks_test = pickle.load(file)

    return np.load(bin_path), np.load(count_path), ks_test

def load_torsions_from_npz(arguments):
    archive = os.path.join(arguments.job, 'histogram', 'torsions.npz')
    return np.load(archive)

def load_torsions_from_dat(arguments):
    all_torsions = {}
    input_path = os.path.join(arguments.job, 'coordinates.dat')
    output_path = os.path.join(arguments.job, 'histogram', 'torsions.npz')
    #progress_bar = '\033[?25l\r[%%d/%s]' % parse_iterations(arguments)
    progress_bar = '\r[%%d/%s]' % parse_iterations(arguments)

    # Read the torsion data off of disk an into standard lists.
    try:
        with open(input_path) as file:
            for line in file:
                try:
                    field, data = line.strip().split(' ', 1)
                except ValueError:
                    continue

                # Use the iteration fields to keep track of progress.
                if field == 'iteration':
                    iteration = int(data)
                    if arguments.limit and iteration > arguments.limit:
                        break
                    else:
                        sys.stdout.write(progress_bar % iteration)
                        sys.stdout.flush()
                        continue

                # Use the phi, psi, omega, and chi fields to extract data.
                elif field in ('phi', 'psi', 'omega', 'chi'):
                    torsions = [float(x) for x in data.split()]
                    for index, value in enumerate(torsions):
                        key = key_from_torsion(field, index)
                        all_torsions.setdefault(key, []).append(value)
    finally:
        #sys.stdout.write('\033[?25h\n')
        sys.stdout.write('\n')

    # Convert the torsion data into numpy arrays.
    for key in all_torsions:
        scratchpad = np.array(all_torsions[key])
        scratchpad %= 360
        scratchpad[scratchpad > 180] -= 360
        all_torsions[key] = scratchpad

    # Write the torsion data to disk.
    np.savez_compressed(output_path, **all_torsions)

    return all_torsions
    
def generate_histograms(arguments, all_torsions):
    import scipy.stats

    bin_path = os.path.join(arguments.job, 'histogram', 'bins.npz')
    count_path = os.path.join(arguments.job, 'histogram', 'counts.npz')
    ks_path = os.path.join(arguments.job, 'histogram', 'ks_test.pkl')

    bins = {}
    counts = {}
    ks_test = {}

    for key in all_torsions:
        torsions = all_torsions[key]

        # Calculate a histogram for each torsion angle.  Each histogram is 
        # characterized by a bins (x-axis) and counts (y-axis) array.  720 bins 
        # are always generated, but the plotting code may display fewer.

        raw_bins = np.linspace(-180, 180, 721)
        counts[key] = np.histogram(torsions, bins=raw_bins)[0]

        # The bins fed into histogram() delimit where each bin starts and 
        # stops.  That's nice, but I need to know where the middle of each bin 
        # is.  These lines convert bin boundaries to bin centers.

        offset = (raw_bins[1] - raw_bins[0]) / 2
        bins[key] = raw_bins[:-1] + offset

        # Calculate a uniformity statistic for each histogram.  A uniform 
        # distribution is expected for linear peptides with the breadth move 
        # enabled.  Use the '-ks' argument to see these results.

        expected = scipy.stats.uniform(loc=-180, scale=360).cdf
        ks_test[key] = scipy.stats.kstest(torsions, expected)[1]

    # Cache these results for later use, then return.

    np.savez_compressed(bin_path, **bins)
    np.savez_compressed(count_path, **counts)

    with open(ks_path, 'w') as file:
        pickle.dump(ks_test, file)

    return bins, counts, ks_test


figure = plt.figure()

def find_title(arguments):
    title_path = os.path.join(arguments.job, 'title.dat')

    if arguments.title:
        return arguments.title

    elif os.path.exists(title_path):
        with open(title_path) as file:
            return file.readline().strip()

    elif arguments.job == '.':
        directories = os.getcwd().split('/')
        return '/'.join(directories[-2:])

    else:
        return arguments.job

def find_y_limits(arguments, axes):
    ylim_path = os.path.join(arguments.job, 'ylim.dat')
    ylims = axes.get_ylim()

    if arguments.ylim:
        ylims = tuple(float(x) for x in arguments.ylim)

    elif os.path.exists(ylim_path):
        with open(ylim_path) as file:
            ylims = tuple(float(x) for x in file.readline().split())

    return ylims[0] - 0.03 * ylims[1], ylims[1]

def find_pivots(arguments):
    path = os.path.join(arguments.job, 'pivots.dat')

    with open(path) as file:
        raw_indices = [int(x) - 1 for x in file.readline().split()]
        pivot_indices = tuple(x - raw_indices[0] for x in raw_indices)

    pivot_keys = []

    for index in pivot_indices:
        pivot_keys.append(key_from_torsion('phi', index))
        pivot_keys.append(key_from_torsion('psi', index))

    return pivot_keys

def bin_histograms(arguments, key, all_bins, all_counts, warning_messages):
    """ Return bins and counts for the given key, where the number of bins is 
    adjusted to best match the value specified by the '--bins' argument.  This 
    is necessary because the same number of bins (720) are always cached, but 
    the user may want a coarser view of the data. """

    fine_bins = all_bins[key]
    fine_counts = all_counts[key]
    fine_size = len(fine_bins)

    # If more bins are requested than are available in the fine arrays, just 
    # return the fine arrays and print a message explaining that higher 
    # resolution data is not available.

    if arguments.bins > fine_size:
        warning_message = "Only using %d bins; no higher resolution available."
        warning_messages.add(warning_message % fine_size)
        return fine_bins, fine_counts

    # The number of bins must divide evenly into the number of stored bins.  If 
    # the requested number of bins is does not evenly divide the number of 
    # stored bins, find the closest number that does.

    best_divisor = None
    best_distance = float('inf')

    for divisor in range(1, fine_size + 1):
        if fine_size % divisor == 0:
            distance = abs(arguments.bins - divisor)

            if distance < best_distance:
                best_divisor = divisor
                best_distance = distance

    divisor = best_divisor
    chunk = fine_size // divisor

    if best_distance != 0:
        warning_message = "Using %d bins instead of %d."
        warning_messages.add(warning_message % (divisor, arguments.bins))

    # Condense the bins and counts to match the requested bin size.

    bins = np.zeros(divisor)
    counts = np.zeros(divisor)

    for index in range(fine_size):
        bins[index//chunk] += fine_bins[index] / chunk
        counts[index//chunk] += fine_counts[index]

    if arguments.normalize:
        counts /= sum(counts)

    return bins, counts

def plot_histograms(arguments, all_bins, all_counts, ks_test, index):
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['axes.titlesize'] = 18
    if os.uname() != 'alanine':
        #mpl.rcParams['font.family'] = 'serif'
        #mpl.rcParams['font.serif'] = 'Liberation Serif'
        mpl.rcParams['font.family'] = 'sans'
        mpl.rcParams['font.serif'] = 'Dejavu Sans'

    # Load the user-specified style.
    num_styles = len(arguments.style)
    module = 'styles.%s' % arguments.style[index % num_styles]
    style = __import__(module, fromlist=['style']).style

    if index == 0: label = "uniform"
    if index == 1: label = "balanced-rama"
    if index == 2: label = "naive-rama"

    # Create a figure which space for an external legend.
    axes = plt.subplot(111)
    box = axes.get_position()

    if arguments.ks_test:
        axes.set_position([box.x0, box.y0, box.width * 0.72, box.height])
    else:
        axes.set_position([box.x0, box.y0, box.width * 0.85, box.height])

    # Plot all the requested histograms.
    warning_messages = set()
    pivots = find_pivots(arguments)

    def torsion_key(key):
        name, index = torsion_from_key(key)
        torsion_order = dict(phi=1, psi=2, omega=3, chi=4)
        return index, torsion_order[name]


    for key in sorted(all_bins, key=torsion_key):
        name, index = torsion_from_key(key)
        uniformity = ks_test[key]
        index = str(index)

        # Get rid of histograms that aren't in the 'include' list.
        if not arguments.include: pass
        elif 'pivot' in arguments.include and key in pivots: pass
        elif 'nonpivot' in arguments.include and key not in pivots: pass
        elif (key not in arguments.include) and \
           (name not in arguments.include) and \
           (index not in arguments.include): continue

        # Get rid of histograms that are in the 'exclude' list.
        if key in arguments.exclude: continue
        if name in arguments.exclude: continue
        if index in arguments.exclude: continue
        if 'pivot' in arguments.exclude and key in pivots: continue
        if 'nonpivot' in arguments.exclude and key not in pivots: continue

        bins, counts = bin_histograms(
                arguments, key, all_bins, all_counts, warning_messages)

        if key in pivots:
            name = name.title()

        if arguments.ks_test:
            attributes = dict(label=r'$\{0}_{1}$ ({2:1.2g})'.format(
                name, index, uniformity))
        else:
            attributes = dict(label=r'$\{0}_{1}$'.format(name, index))

        #attributes['label'] += ' ' + label

        attributes.update(style.get(key, {}))

        axes.plot(bins, counts, **attributes)

    for message in warning_messages:
        print "Warning:", message

    # Format the graph.
    title = find_title(arguments)
    ylim = find_y_limits(arguments, axes)

    axes.set_title(title)
    axes.set_xlabel('Torsion DOF')
    axes.set_xticks(range(-180, 181, 60))
    axes.set_xlim(-180, 180)
    axes.set_ylim(*ylim)
    axes.legend(loc='upper left', bbox_to_anchor=(1.05, 1))

def plot_score_function(arguments):
    from graphics import tango

    if not arguments.rama:
        return

    data = {}

    with open('data/rama_score_function.dat') as file:
        for line in file:
            phi, psi, score = line.split()
            phi, psi, score = int(phi), int(psi), float(score)
            data.setdefault(phi, {})[psi] = score

    phis = np.sort(data.keys())
    psis = np.sort(data.values()[0].keys())
    scores = np.zeros((len(phis), len(psis)))

    for i, phi in enumerate(phis):
        for j, psi in enumerate(psis):
            scores[i,j] = data[phi][psi]

    phi_scores = np.exp(-scores).sum(1)
    psi_scores = np.exp(-scores).sum(0)

    phi_scores /= np.sum(phi_scores)
    psi_scores /= np.sum(psi_scores)

    plt.plot(phis, phi_scores, '-', color=tango.grey[1], lw=10.0)
    plt.plot(psis, psi_scores, '-', color=tango.grey[1], lw=10.0)


if __name__ == '__main__':
    arguments = parse_arguments()

    plot_score_function(arguments)

    for index, job in enumerate(arguments.jobs):
        arguments.job = job
        bins, counts, ks_test = load_histograms(arguments)
        plot_histograms(arguments, bins, counts, ks_test, index)

    # Output the graph in the desired format.
    if arguments.output:
        plt.savefig(arguments.output, dpi=arguments.dpi)
    if not arguments.no_gui and not os.fork():
        plt.show()


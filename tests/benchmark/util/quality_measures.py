#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   benchmark.py
## @brief  Helper functions for quality measures
## @author Sergey Lyskov

import os, os.path, sys, imp, shutil, json, math
import numpy as np

#=======================================
def check_all_values_above_cutoff( col, cutoff, tag, filehandle ):

	out = "All " + tag + "s > cutoff"
	filehandle.write( out + " " + str(cutoff) + "\t" )

	if all( i >= cutoff for i in col ):
		value = True
	else:
		value = False

	filehandle.write( str( value ) + "\n" )
	return {out : value}

#=======================================
def check_all_values_below_cutoff( col, cutoff, tag, filehandle ):

	return check_xpercent_values_below_cutoff( col, cutoff, tag, filehandle, 100 )

# =======================================
def kabsch_align(r1_coords: np.ndarray, r2_coords: np.ndarray):
    """Perform kabsch align on two sets of coordinates
    Notes:
        We align two sets of numpy coordinates (indicies must match, ie you must want
        idx 0 of 'r1_coords1' and 'r2_coords' to align together).
    """
    E0 = np.sum(np.sum(r1_coords * r1_coords, axis=0), axis=0) + np.sum(np.sum(r2_coords * r2_coords, axis=0), axis=0)
    V, S, Wt = np.linalg.svd(np.dot(np.transpose(r2_coords), r1_coords))
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))
    if np.isclose(-1.0, reflect):
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    rotation_matrix = np.dot(V, Wt)
    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / len(r1_coords)))
    return rotation_matrix, RMSD


# =======================================
def check_all_values_below_cutoff(col, cutoff, tag, filehandle):

    return check_xpercent_values_below_cutoff(col, cutoff, tag, filehandle, 100)


# =======================================
def check_xpercent_values_below_cutoff(col, cutoff, tag, filehandle, percentage):

    out = "All " + tag + "s < cutoff"
    filehandle.write(out + " " + str(cutoff) + "\t")

    # sort the values from smallest to largest, then take the first x records
    col = sorted(col)
    partial_col = col[: int(math.ceil(percentage * len(col) / 100.0))]

    if all(i <= cutoff for i in partial_col):
        value = True
    else:
        value = False

    filehandle.write(str(value) + "\n")
    return {out: value}


# =======================================
def calc_Conway_discrim_score(rmss, scores, offsets=[0, 0.5, 1, 1.5, 2, 4]):
    # calculate funnel discrimination metric from Conway 2014
    # "Relaxation of backbone bond geometry improves protein
    # energy landscape modeling"
    # with a slight twist, rms cuts are based on min rms + an offset

    scores = norm_Conway(scores)
    rmss = np.array(rmss)  # must be np.array for indexing

    D = 0.0  # discrimination score
    cuts = [min(rmss) + x for x in offsets]
    for r in cuts:
        # lowest score below rms threshold is set to 0 to avoid NAs
        below = 0
        if len(scores[rmss <= r]) > 0:
            below = min(scores[rmss <= r])
        above = below  # get lowest score above rms threshold, default to below so diff is 0 if NAs
        if len(scores[rmss > r]) > 0:
            above = min(scores[rmss > r])
        D += below - above  # difference in scores across thresholds is the disrcim

    return round(D, 3)


# =======================================
def norm_Conway(scores):
    # from Conway 2014 "Relaxation of backbone bond geometry
    # improves protein energy landscape modeling"
    # normalize scores 95th percent = 1, 5th percent = 0

    s_scores = sorted(scores)  # sort low to high

    high_score = s_scores[int(0.95 * len(s_scores))]
    low_score = s_scores[int(0.05 * len(s_scores))]

    scores = (np.array(scores) - low_score) / (high_score - low_score)

    return scores


# =======================================
def check_rmsd_of_topscoring(rmsd_col_sorted, cutoff, filehandle):

    out = "rmsd of topscoring model < " + str(cutoff)
    filehandle.write(out + "\t")

    if rmsd_col_sorted[0] <= cutoff:
        value = True
    else:
        value = False

    filehandle.write(str(value) + "\n")
    return {out: value}


# =======================================
def check_range(col, tag, filehandle):

    filehandle.write(
        tag
        + "\tmin, max, avg, std:"
        + "% 12.3f % 12.3f % 12.3f % 12.3f\n" % (min(col), max(col), np.mean(col), np.std(col))
    )
    value = {
        "min": round(min(col), 4),
        "max": round(max(col), 4),
        "avg": round(np.mean(col), 4),
        "std": round(np.std(col), 4),
    }
    return {tag: value}


# =======================================
def check_for_2d_top_values(
    xcol, ycol, xtag, ytag, xcutoff, ycutoff, filehandle, xminimize=True, yminimize=True, xabsolute=True, yabsolute=True
):
    """
    Check if there are decoys that are in both the top xcutoff% and top ycutoff%.
    Returns a tuple, containing (number of passes, x_absolute_cutoff, y_absolute_cutoff).
    xcol, ycol = Arrays of scoreterm values to be checked
    xtag, ytag = Used to genereate the out tag.  Should usually be the scoreterm names being checked (but is arbitrary)
    xcutoff, ycutoff = The absolute or percentile cutoff for xcol and ycol, respectively.  If percentile, should be between 0 and 100.
    filehandle = The result file to write to.
    xminimize, yminimize = By default (True), the "top" values are those less than the cutoff.  If false, the "top" values are those greater than the cutoff.
    xabsolute, yabsolute = Should the threshold be using an absolute cutoff (True) or a percentile cutoff (False).
    """
    assert len(xcol) == len(ycol)

    out = "Num decoys with " + xtag + " and " + ytag + " better than cutoffs"
    filehandle.write(out + " " + str(xcutoff) + " and " + str(ycutoff) + ":\t")

    # Set the cutoff values in x and y
    x_val = xcutoff
    y_val = ycutoff
    # Get the percentile values for xcol and ycol if we ask for it
    if not xabsolute:
        x_val = np.percentile(xcol, xcutoff)
    if not yabsolute:
        y_val = np.percentile(ycol, ycutoff)

    # Check each value in xcol to see if it is greater/less than xperc_val
    # Then check if the corresponding value in ycol is greater/less than yperc_val
    # If both passing conditions are met, append True to passed.  If either condition fails, append False to passed.
    passed = []
    for idx in range(len(xcol)):
        if xminimize:
            if xcol[idx] <= x_val:
                if yminimize:
                    if ycol[idx] <= y_val:
                        passed.append(True)
                    else:
                        passed.append(False)
                else:  # ymaximize
                    if ycol[idx] >= y_val:
                        passed.append(True)
                    else:
                        passed.append(False)
            else:
                passed.append(False)
        else:  # xmaximize
            if xcol[idx] >= x_val:
                if yminimize:
                    if ycol[idx] <= y_val:
                        passed.append(True)
                    else:
                        passed.append(False)
                else:  # ymaximize
                    if ycol[idx] >= y_val:
                        passed.append(True)
                    else:
                        passed.append(False)
            else:
                passed.append(False)
    assert len(passed) == len(xcol)
    assert len(passed) == len(ycol)

    # The number of passes is the sum of passed (true = 1, false = 0)
    npasses = sum(passed)
    filehandle.write(str(npasses) + " (" + str(x_val) + ", " + str(y_val) + ")\n")
    return {out: (npasses, x_val, y_val)}

#=======================================
# @brief Given a list of scores, a matching list of rmsds, and values for lambda and kbt,
# compute the PNear value (a measure of fold propensity, with 0 meaning that the sequence does
# not spend any time in the native state, and 1 meaning that it spends all of its time in the
# native state).
# @details Lambda is a value in Angstroms indicating the breadth of the Gaussian used to define
# "native-like-ness".  The bigger the value, the more permissive the calculation is to structures
# that deviate from native.  Typical values for peptides range from 1.5 to 2.0, and for proteins
# from 2.0 to perhaps 4.0.  The value of kbt, in energy units, determines how large an energy gap
# must be in order for a sequence to be said to favour the native state.  The default value, 0.62,
# should correspond to physiological temperature for ref2015 or any other scorefunction with units
# of kcal/mol.
# @note Unlike the Conway discrimination score, the PNear calculation uses no hard cutoffs.  This is
# advantageous for repeated testing: if the scatter of points on the RMSD plot changes very slightly
# from run to run, the PNear value will only change by a small amount, whereas any metric dependent
# on hard cutoffs could change by a large amount if a low-energy point crosses an RMSD threshold.
# @author Vikram K. Mulligan (vmulligan@flatironinstitute.org
def calculate_pnear( scores, rmsds, lambda_val=1.5, kbt=0.62 ) :
    nscores = len(scores)
    assert nscores == len(rmsds), "Error in calculate_pnear(): The scores and rmsds lists must be of the same length."
    assert nscores > 0, "Error in calculate_pnear(): At least one score/rmsd pair must be provided."
    assert kbt > 1e-15, "Error in calculate_pnear(): kbt must be greater than zero!"
    assert lambda_val > 1e-15, "Error in calculate_pnear(): lambda must be greater than zero!"
    minscore = min( scores )
    weighted_sum = 0.0
    Z = 0.0
    lambdasq = lambda_val * lambda_val
    for i in range( nscores ) :
        val1 = exp( -( rmsds[i] * rmsds[i] ) / lambdasq )
        val2 = exp( -( scores[i] - minscore ) / kbt )
        weighted_sum += val1*val2
        Z += val2
    assert Z > 1e-15, "Math error in calculate_pnear()!  This shouldn't happen."
    return weighted_sum/Z

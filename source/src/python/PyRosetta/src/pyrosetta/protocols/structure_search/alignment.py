import collections

import numpy

import pyrosetta.numeric.alignment.rmsd_calc as rmsd_calc


def length_blocks(jagged_array):
    blocks = collections.defaultdict(list)

    for i, l in enumerate(map(len, jagged_array)):
        blocks[l].append(i)

    return dict(blocks)


def length_grouped_pairwise_metric(samples, pairwise_metric):
    if not isinstance(samples, numpy.ndarray):
        samples = numpy.array(samples)

    if samples.ndim == 2:
        return pairwise_metric(samples)
    else:
        assert samples.ndim == 1

        result = numpy.full((len(samples), len(samples)), numpy.nan, float)

        for l, ind in length_blocks(samples).items():
            sub = length_grouped_pairwise_metric(list(samples[ind]), pairwise_metric)
            for i in range(len(ind)):
                result[ind[i]][ind] = sub[i]

        return result


def length_grouped_broadcast_metric(samp_a, samp_b, broadcast_metric):
    if not isinstance(samp_a, numpy.ndarray):
        samp_a = numpy.array(samp_a)
    if not isinstance(samp_b, numpy.ndarray):
        samp_b = numpy.array(samp_b)

    if samp_a.ndim == samp_b.ndim == 2 and samp_a.shape[-1] == samp_b.shape[-1]:
        return broadcast_metric(samp_a, samp_b)

    result = numpy.full((len(samp_a), len(samp_b)), numpy.nan, float)

    b_inds = length_blocks(samp_b)
    b_blocks = {l: numpy.array(list(samp_b[inds])) for l, inds in b_inds.items()}

    for a in range(len(samp_a)):
        al = len(samp_a[a])

        result[a][b_inds[al]] = broadcast_metric(
                numpy.array(samp_a[a]).reshape((1, -1)), b_blocks[al])

        return result


def fragment_unaligned_pairwise_rmsd(fragment_residues):
    """Calculate pairwise rmsd over given, potentially uniqual length, residue arrays."""
    def unaligned_rmsd(res):
        return rmsd_calc.structured_array_unaligned_broadcast_rmsd(res["orient"], res["orient"])

    return length_grouped_pairwise_metric(fragment_residues, unaligned_rmsd)


def fragment_pairwise_rmsd(fragment_residues):
    """Calculate pairwise rmsd over given, potentially uniqual length, residue arrays."""
    def aligned_rmsd(res):
        return rmsd_calc.structured_array_broadcast_rmsd(res["orient"], res["orient"])

    return length_grouped_pairwise_metric(fragment_residues, aligned_rmsd)


def fragment_unaligned_broadcast_rmsd(fragment_residues_a, fragment_residues_b):
    """Calculate unaligned rmsd over given, potentially uniqual length, residue arrays."""
    def unaligned_rmsd(ra, rb):
        return rmsd_calc.structured_array_unaligned_broadcast_rmsd(ra["orient"], rb["orient"])

    return length_grouped_broadcast_metric(fragment_residues_a, fragment_residues_b, unaligned_rmsd)


def fragment_broadcast_rmsd(fragment_residues_a, fragment_residues_b):
    """Calculate aligned rmsd over given, potentially uniqual length, residue arrays."""
    def aligned_rmsd(ra, rb):
        return rmsd_calc.structured_array_broadcast_rmsd(ra["orient"], rb["orient"])

    return length_grouped_broadcast_metric(
            fragment_residues_a, fragment_residues_b, aligned_rmsd)

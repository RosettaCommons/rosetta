"""Utility functions related to rosetta energies."""

try:
    import numpy
except ImportError:
    pass

__all__ = [
    "residue_pair_energies_array",
    "residue_onebody_energies_array",
    "residue_total_energies_array",
    "total_energies_array",
    "nonzero_weights",
    "nonzero_weights_dtype",
    "energies_total_score",
]

from pyrosetta.bindings.utility import bind_method
import pyrosetta.rosetta.core.scoring


def _residue_selection_to_1ind(selection, total_count):
    if selection is None:
        selection = numpy.arange(total_count) + 1

    selection = numpy.array(selection)
    if selection.dtype == bool:
        selection = numpy.flatnonzero(selection) + 1
    else:
        selection = selection.astype(int)

    assert numpy.alltrue(selection > 0), "Selection contained invalid indicies."
    assert numpy.alltrue(selection <= total_count), "Selection contained invalid indicies."

    return selection


@bind_method(pyrosetta.rosetta.core.scoring.Energies)
def residue_pair_energies_array(energies, from_residue_selection=None, to_residue_selection=None):
    """Generate pair energy table from the given energies object.
    returns: energy_table - 2d energy type structured array.
                shape - (energies.size, energies.size)
                dtype - [(<type>, float)] for every nonzero weight energy type.
    """

    energy_graph = energies.energy_graph()

    from_residue_selection = _residue_selection_to_1ind(
        from_residue_selection, energy_graph.num_nodes())
    to_residue_selection = _residue_selection_to_1ind(
        to_residue_selection, energy_graph.num_nodes())

    energy_types = list(energy_graph.active_2b_score_types())
    energy_table = numpy.zeros(
        (len(from_residue_selection), len(to_residue_selection), len(energy_types)))

    for i in range(len(from_residue_selection)):
        for j in range(len(to_residue_selection)):
            edge = energy_graph.find_edge(from_residue_selection[i], to_residue_selection[j])
            if edge:
                for t in range(len(energy_types)):
                    energy_table[i, j, t] = edge[energy_types[t]]

    pair_term_table = energy_table.view(
        dtype=[(str(etype).split(".")[-1], float) for etype in energy_types]
    ).reshape(energy_table.shape[:-1])

    pair_table = numpy.zeros_like(
        pair_term_table,
        dtype=[("total_score", float)] + pair_term_table.dtype.descr
    )

    for n in pair_term_table.dtype.names:
        pair_table[n] = pair_term_table[n]

    pair_table["total_score"] = energies_total_score(pair_term_table, nonzero_weights(energies))

    return pair_table


@bind_method(pyrosetta.rosetta.core.scoring.Energies)
def residue_onebody_energies_array(energies, residue_selection=None, out=None):
    """Gets table of energy terms with non-zero weight on a per-residue basis.
    returns:
        structure_array of shape n_residue, with per-score-term entries
    """

    from pyrosetta.rosetta.core.scoring import ScoreType

    residue_selection = _residue_selection_to_1ind(residue_selection, energies.size())

    if out is None:
        out = numpy.empty(len(residue_selection), nonzero_weights_dtype(energies.weights()))
    else:
        out = out.ravel()
        assert out.shape == (len(residue_selection),)

    scoretypes = [
        (name, ScoreType.__dict__[name])
        for name in out.dtype.names if name != "total_score"
    ]

    for i in range(len(residue_selection)):
        residue_totals = energies.onebody_energies(residue_selection[i])
        for n, t in scoretypes:
            out[i][n] = residue_totals[t]

        out[i]["total_score"] = (residue_totals * energies.weights()).sum()

    return out


@bind_method(pyrosetta.rosetta.core.scoring.Energies)
def residue_total_energies_array(energies, residue_selection=None, out=None):
    """Gets table of energy terms with non-zero weight on a per-residue basis.
    returns:
        structure_array of shape n_residue, with per-score-term entries
    """

    from pyrosetta.rosetta.core.scoring import ScoreType

    residue_selection = _residue_selection_to_1ind(residue_selection, energies.size())

    if out is None:
        out = numpy.empty(len(residue_selection), nonzero_weights_dtype(energies.weights()))
    else:
        out = out.ravel()
        assert out.shape == (len(residue_selection),)

    scoretypes = [(name, ScoreType.__dict__[name]) for name in out.dtype.names]

    for i in range(len(residue_selection)):
        residue_totals = energies.residue_total_energies(residue_selection[i])
        for n, t in scoretypes:
            out[i][n] = residue_totals[t]

    return out


@bind_method(pyrosetta.rosetta.core.scoring.Energies)
def total_energies_array(energies, out=None):
    """Get total structured dtype with non-zero energies."""
    from pyrosetta.rosetta.core.scoring import ScoreType

    if out is None:
        out = numpy.empty(1, nonzero_weights_dtype(energies.weights()))
    else:
        out = out.ravel()
        assert out.shape == (1,)

    scoretypes = [(name, ScoreType.__dict__[name]) for name in out.dtype.names]

    total_energies = energies.total_energies()

    for n, t in scoretypes:
        out[n] = total_energies[t]

    return out

@bind_method(pyrosetta.rosetta.core.scoring.Energies)
def active_total_energies(energies):
    if energies.weights().sum() == 0:
        return {}

    total_array = energies.total_energies_array()
    return { n : float(total_array[n]) for n in total_array.dtype.names }

@bind_method(pyrosetta.rosetta.core.scoring.Energies)
def nonzero_weights(energies, out=None):
    """Gets energy terms weights in the given energies object.
    returns:
        { score_term : score_weight }
    """

    if out is None:
        out = numpy.empty(1, nonzero_weights_dtype(energies.weights()))
    else:
        out = out.ravel()
        assert out.shape == (1,)

    from pyrosetta.rosetta.core.scoring import ScoreType

    scoretypes = [(name, ScoreType.__dict__[name]) for name in out.dtype.names]

    weights = energies.weights()

    for n, t in scoretypes:
        out[n] = weights[t]

    return out


def nonzero_weights_dtype(weights):
    from pyrosetta.rosetta.core.scoring import ScoreType
    return numpy.dtype([
        (name, float) for (name, st) in ScoreType.__dict__.items()
        if not name.startswith("_") and (weights.get(st) != 0 or name == "total_score")
    ])


def energies_total_score(energy_table, energy_weights):
    """Convert energy term struct array into total score array via given weights.
    energies - Struct array with energy term fields.
    energy_weights - dict of energy weights.
    returns - total_score float array of energies.shape.
    """

    total = numpy.zeros_like(energy_table, dtype=float)

    for eterm in energy_table.dtype.names:
        total += energy_table[eterm] * energy_weights[eterm]

    return total

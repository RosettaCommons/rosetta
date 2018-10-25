import numpy
import collections

from .search import (
    StructureDatabase,
    StructurePairQuery,
    PairQueryExecutor,
    StructureSingleQuery,
    SingleQueryExecutor,
)

from pyrosetta.utility.array import structured_array_to_basic


class StructureSearchManager(object):
    """Docstring for StructureSearchManager. """

    def __init__(self, source_residues):
        """TODO: to be defined1.

        Args:
            source_residues (TODO): TODO


        """
        self.source_residues = source_residues
        self.target_atoms = source_residues["orient"].dtype.names

        source_coords = structured_array_to_basic(self.source_residues["orient"])
        structure_endpoints = numpy.flatnonzero(self.source_residues["structure_id"][:-1] != self.source_residues["structure_id"][1:])
        chain_endpoints = numpy.flatnonzero(self.source_residues["chain_ending"])

        self.db = StructureDatabase()
        self.db.initialize(source_coords, structure_endpoints, chain_endpoints)

    def pair_query(self, query_components, query_tolerance, primary_range=None):
        if primary_range:
            if isinstance(primary_range, collections.Iterable):
                min_primary_distance, max_primary_distance = primary_range
            elif primary_range:
                min_primary_distance, max_primary_distance = 0, primary_range

            q = StructurePairQuery(
                structured_array_to_basic(query_components[0][list(self.target_atoms)]),
                structured_array_to_basic(query_components[1][list(self.target_atoms)]),
                query_tolerance,
                min_primary_distance, max_primary_distance
            )
        else:
            q = StructurePairQuery(
                structured_array_to_basic(query_components[0][list(self.target_atoms)]),
                structured_array_to_basic(query_components[1][list(self.target_atoms)]),
                query_tolerance,
            )

        e = PairQueryExecutor(q)
        e.execute(self.db)
        return e.query_results

    def single_query(self, query_component, query_tolerance):
        q = StructureSingleQuery(
            structured_array_to_basic(query_component[list(self.target_atoms)]),
            query_tolerance
        )

        e = SingleQueryExecutor(q)
        e.execute(self.db)
        return e.query_results

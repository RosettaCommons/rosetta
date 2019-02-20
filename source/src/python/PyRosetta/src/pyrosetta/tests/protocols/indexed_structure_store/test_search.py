# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import numpy
import pyrosetta.protocols.structure_search as search
import pyrosetta.numeric.alignment.rmsd_calc as rmsd_calc
import unittest

from collections import namedtuple
from pyrosetta.tests.protocols.indexed_structure_store.data import source_structure_residues
from pyrosetta.utility.array import atom_array_to_coordinates


def fragment_start_indicies_from_residue_array(fragment_length, source_residues):
    chain_endings = source_residues["chain_ending"].copy()
    chain_endings[-1] = True
    chain_endings[:-1][source_residues[:-1]["structure_id"] != source_residues[1:]["structure_id"]] = True
    chain_endings = numpy.flatnonzero(chain_endings)

    fragment_starts = numpy.ones_like(source_residues, dtype=bool)
    for i in range(1, fragment_length):
        fragment_starts[chain_endings - (i - 1)] = False
        fragment_starts[-i] = False

    return numpy.flatnonzero(fragment_starts)


def fragments_from_source_residues(fragment_length, source_residues):
    """Extract fragment entries from residue array.

    source_residues - Source residue data, include id, resn, and 'orient' entry containing atomic coordinates.

    returns - (fragment_start_indices, fragment_data)
    """

    fragment_start_indicies = fragment_start_indicies_from_residue_array(fragment_length, source_residues)

    # Pick source residues from residue store
    fragment_residues = numpy.empty(fragment_start_indicies.shape + (fragment_length, ), source_residues.dtype)
    for i in range(fragment_length):
        fragment_residues[:, i] = source_residues[fragment_start_indicies + i]

    return (fragment_start_indicies, fragment_residues)


class TestStructureSearch(unittest.TestCase):

    def gen_test_set(self):
        # Add duplicate structure to test structure isolation
        dup_res = source_structure_residues[source_structure_residues["structure_id"] == 0].copy()
        dup_res["structure_id"] = 1663
        return numpy.concatenate((
            source_structure_residues[source_structure_residues["structure_id"] < 3], dup_res))

    @unittest.skipIf(source_structure_residues is None, "Missing source data")
    def test_single_search(self):
        test_res = self.gen_test_set()

        test_length = 9
        fstarts, fdata = fragments_from_source_residues(test_length, test_res)
        fcoords = atom_array_to_coordinates(fdata["orient"])

        q_i = 4
        query_threshold = .4
        pd = rmsd_calc.coordinate_array_rmsd(fcoords[q_i], fcoords)
        query_matches = pd < query_threshold

        query_coords = test_res[fstarts[q_i]:fstarts[q_i] + test_length]["orient"]

        # Assert that test data produced an expected result
        assert numpy.count_nonzero(query_matches) > 10

        sm = search.StructureSearchManager(test_res)
        qr = sm.single_query(query_coords, query_threshold)

        self.assertEqual(len(qr), numpy.count_nonzero(query_matches))
        self.assertEqual(set(qr["fragment_start"]), set(fstarts[query_matches]))

    @unittest.skipIf(source_structure_residues is None, "Missing source data")
    def test_pair_search(self):
        test_res = self.gen_test_set()

        test_length = 2
        fstarts, fdata = fragments_from_source_residues(test_length, test_res)
        fcords = atom_array_to_coordinates(fdata["orient"])

        ffrom, fto = numpy.nonzero(fdata[:, 0]["structure_id"].reshape((-1, 1)) == fdata[:, 0]["structure_id"].reshape((1, -1)))

        pcoords = numpy.concatenate((fcords[ffrom], fcords[fto]), axis=1)
        assert pcoords.shape == (len(ffrom), 16, 3)

        # Query index is index of query fragment pair in array above
        q_i = 4
        # Query fragment indicies in test_res
        qs1, qs2 = fstarts[ffrom[q_i]], fstarts[fto[q_i]]
        query_threshold = .75
        query_coords = test_res[qs1:qs1 + 2]["orient"], test_res[qs2:qs2 + 2]["orient"]
        pd = rmsd_calc.coordinate_array_rmsd(pcoords[q_i], pcoords)
        query_matches = pd < query_threshold

        # Assert that test data produced an expected result
        assert numpy.count_nonzero(query_matches) > 10
        assert all(fstarts[fto][query_matches] - fstarts[ffrom][query_matches] == qs2 - qs1)

        sm = search.StructureSearchManager(test_res)
        qr = sm.pair_query(query_coords, query_threshold)

        self.assertEqual(len(qr), numpy.count_nonzero(query_matches))
        self.assertEqual(
            set(zip(qr["fragment_a_start"], qr["fragment_b_start"])),
            set(zip(fstarts[ffrom][query_matches], fstarts[fto][query_matches]))
        )

    @unittest.skipIf(source_structure_residues is None, "Missing source data")
    def test_rep_search(self):
        src_res = source_structure_residues[source_structure_residues["structure_id"] == 0].copy()

        start_fragment = src_res[4:6]
        int_res = [src_res[6:7]]
        end_fragment = src_res[7:9]

        structures = []
        for il in range(0, 10):
            s = numpy.concatenate([start_fragment] + int_res * il + [end_fragment])
            s["structure_id"] = il
            structures.append(s)

        sub_results = []
        for i, s in enumerate(structures):
            search_manager = search.StructureSearchManager(s)
            search_results = search_manager.pair_query((start_fragment["orient"], end_fragment["orient"]), .5)

            self.assertEqual(len(search_results), 1)
            self.assertEqual(search_results[0]["fragment_a_start"], 0)
            self.assertEqual(search_results[0]["fragment_b_start"], len(start_fragment) + i)
            sub_results.append(search_results)

        test_res = numpy.concatenate(structures)
        tdb = search.StructureSearchManager(test_res)
        tr = tdb.pair_query((start_fragment["orient"], end_fragment["orient"]), .5)

        self.assertEqual(len(tr), len(sub_results))

        soffset = 0
        for i, s in enumerate(structures):
            # Expected inter-fragment offset
            el = len(start_fragment) + i
            sr = tdb.pair_query((start_fragment["orient"], end_fragment["orient"]), .5, (el, el))
            self.assertEqual(len(sr), 1)
            self.assertEqual(sr[0]["fragment_a_start"], soffset)
            self.assertEqual(sr[0]["fragment_b_start"], soffset + len(start_fragment) + i)

            soffset += len(s)

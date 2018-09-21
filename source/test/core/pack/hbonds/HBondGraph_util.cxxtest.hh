// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/hbonds/HBondGraph_util.cxxtest.hh
/// @brief  test suite for HBondGraph_util.hh
/// @author Brian Coventry (bcov@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

#include <platform/types.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <utility/graph/Graph.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>
#include <core/pack/hbonds/HBondGraph_util.hh>


#include <core/types.hh>

#include <test/UTracer.hh>
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;
using namespace core::scoring::hbonds::graph;
using namespace core::pack::hbonds;


static basic::Tracer TR("core.pack.hbonds.HBondGraph_util.cxxtest");

class HBondGraphUtilTests : public CxxTest::TestSuite
{
public:
	void setUp()
	{
		core_init();
	}

	void atom_has_match(
		utility::vector1< AtomInfo > const & atoms,
		utility::vector1< bool > & atom_used,
		conformation::Residue const & res,
		Size at
	) {
		bool matched = false;
		for ( Size iinfo = 1; iinfo <= atoms.size(); iinfo++ ) {
			AtomInfo const & info = atoms.at(iinfo);
			if ( info.local_atom_id() == at ) {
				TS_ASSERT( ! matched );
				matched = true;
				atom_used.at(iinfo) = true;
				TS_ASSERT( info.xyz().distance( res.xyz( at ) ) < 0.1f );
			}
		}
		TS_ASSERT( matched );
	}

	void test_graph_completeness(){
		pose::Pose pose;
		import_pose::pose_from_file( pose, "core/pack/dunbrack/1UBQ_repack.pdb" , core::import_pose::PDB_file);

		// We only need hbonds for this test

		scoring::ScoreFunctionOP sfxn_sc ( new scoring::ScoreFunction() );
		sfxn_sc->set_weight( scoring::hbond_bb_sc, 1.0 );
		sfxn_sc->set_weight( scoring::hbond_sc, 1.0 );
		sfxn_sc->score( pose );

		scoring::ScoreFunctionOP sfxn_bb ( new scoring::ScoreFunction() );
		sfxn_bb->set_weight( scoring::hbond_sr_bb, 1.0 );
		sfxn_bb->set_weight( scoring::hbond_lr_bb, 1.0 );
		sfxn_bb->score( pose );


		// Use hbond_graph_from_partial_rotsets to fill the hb_graph

		pack::rotamer_set::RotamerSetsOP blank_rotsets( new pack::rotamer_set::RotamerSets() );
		pack::task::PackerTaskOP task = pack::task::TaskFactory::create_packer_task( pose );
		task->temporarily_fix_everything();
		blank_rotsets->set_task( task );

		pack::rotamer_set::RotamerSetsOP rotsets = nullptr;
		utility::vector1<bool> not_used;

		scoring::hbonds::graph::AtomLevelHBondGraphOP hb_graph;
		hb_graph = core::pack::hbonds::hbond_graph_from_partial_rotsets(
			pose, *blank_rotsets, sfxn_sc, sfxn_bb, rotsets, not_used, -0.001f );


		// The hb_graph should contain all polar atoms and all hbonds. Lets check

		// First check for all polar atoms
		TS_ASSERT( hb_graph->num_nodes() == pose.size() );

		for ( Size resnum = 1; resnum <= pose.size(); resnum++ ) {
			AtomLevelHBondNode * node = hb_graph->get_node( resnum );
			conformation::Residue const & res = pose.residue(resnum);

			TS_ASSERT( node->moltenres() == rotsets->resid_2_moltenres( resnum ) );

			utility::vector1< AtomInfo > const & atoms = node->polar_sc_atoms_not_satisfied_by_background();
			utility::vector1< bool > atom_used( atoms.size(), false );

			// Make sure all polar heavy atoms appear exactly once
			for ( Size at = 1; at <= res.nheavyatoms(); at++ ) {

				// Use BuriedUnsat style check for cross-redundancy
				if ( res.atom_type( at ).is_acceptor() || res.atom_type( at ).is_donor() ) {
					// have to check for backbone N on proline ( this can go away once Npro no longer incorrectly typed as DONOR
					if ( res.name1() == 'P' && res.atom_type( at ).atom_type_name() == "Npro" ) continue;

					atom_has_match( atoms, atom_used, res, at);

				}
			}

			// Make sure all polar hydrogens appear exactly once
			for ( Size at = res.nheavyatoms() + 1; at <= res.natoms(); at++ ) {
				if ( res.atom_is_polar_hydrogen( at ) ) {

					atom_has_match( atoms, atom_used, res, at);
				}
			}

			// Make sure all hb_graph atoms were actually found
			for ( Size i = 1; i <= atom_used.size(); i++ ) {
				TS_ASSERT( atom_used[i] );
			}
		}

		// Next check for all hbonds
		scoring::hbonds::HBondSet hbset;
		scoring::hbonds::fill_hbond_set( pose, false, hbset );
		utility::vector1< scoring::hbonds::HBondOP >  const & hbonds = hbset.hbonds();
		utility::vector1<bool> hbond_used( hbonds.size(), false );


		// Upper triangle for-loop
		for ( Size ii = 1; ii <= pose.size(); ii++ ) {
			for ( Size jj = 1; jj <= pose.size(); jj++ ) {
				if ( ii >= jj ) continue;

				AtomLevelHBondEdge const * edge = hb_graph->find_edge( ii, jj );
				if ( ! edge ) continue;
				Size first_id = edge->get_first_node_ind();
				Size second_id = edge->get_second_node_ind();

				for ( HBondInfo const & info : edge->hbonds() ) {

					Size donor_res = info.first_node_is_donor() ? first_id : second_id;
					Size acceptor_res = info.first_node_is_donor() ? second_id : first_id;

					bool matched = false;
					for ( Size ihb = 1; ihb <= hbonds.size(); ihb++ ) {
						scoring::hbonds::HBondOP const & hb = hbonds.at(ihb);
						if ( hb->acc_res() != acceptor_res || hb->don_res() != donor_res ) continue;

						if ( info.local_atom_id_A() != hb->acc_atm() ) continue;
						if ( info.local_atom_id_H() != hb->don_hatm() ) continue;

						TS_ASSERT( !matched );
						matched = true;
						hbond_used.at(ihb) = true;
					}
					TS_ASSERT( matched );
				}
			}
		}

		// Make sure all hbset hbonds were actually found
		for ( Size i = 1; i <= hbond_used.size(); i++ ) {
			TS_ASSERT( hbond_used[i] );
		}


	}


};

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/vdwaals/VDW_Energy.cxxtest.hh
/// @brief  Unit tests for the centroid van der Waals term
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/types.hh>

// Unit headers
#include <core/scoring/vdwaals/VDW_Energy.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/util/SwitchResidueTypeSet.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

#include <core/conformation/AbstractRotamerTrie.hh>
#include <core/scoring/trie/RotamerTrieBase.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.etable.EtableEnergy.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// @name EtableEnergyTest
/// @brief: Test the functionality of the EtableEnergy class
///////////////////////////////////////////////////////////////////////////
class VDW_EnergyTests : public CxxTest::TestSuite {

public:

	void setUp() {
		/// The analytic version operates just fine in the presence of the table-based etable, but
		/// the table-based version won't work unless the analytic_etable_evaluation flag is false
		/// at the time the FA_STANDARD Etable is constructed.
		core_init();
	}

	void tearDown() {}


	void test_vdw_trie_vs_trie() {
		using namespace core::graph;
		using namespace core::pose;
		using namespace core::scoring::vdwaals;
		using namespace core::scoring::methods;
		using namespace core::pack;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;

		Pose pose = create_trpcage_ideal_pose();
		core::util::switch_to_residue_type_set( pose, "centroid_rot" );
		ScoreFunction sfxn;
		sfxn.set_weight( vdw, 0.5 );
		sfxn( pose );

		EnergyMethodOptions options; // default is fine
		VDW_Energy vdw_energy( options );

		PackerTaskOP task = TaskFactory::create_packer_task( pose );
		for ( Size ii = 1; ii <= 7; ++ii ) task->nonconst_residue_task(ii).prevent_repacking();
		for ( Size ii = 12; ii <= pose.total_residue(); ++ii ) task->nonconst_residue_task(ii).prevent_repacking();
		for ( Size ii = 8; ii <= 11; ++ii ) task->nonconst_residue_task( ii ).or_ex1( true );
		for ( Size ii = 8; ii <= 11; ++ii ) task->nonconst_residue_task( ii ).or_ex2( true );

		GraphOP packer_neighbor_graph = create_packer_graph( pose, sfxn, task );
		RotamerSets rotsets; rotsets.set_task( task );
		rotsets.build_rotamers( pose, sfxn, packer_neighbor_graph );

		// Create tries for each of the three rotamer sets
		for ( Size ii = 8; ii <= 11; ++ii ) vdw_energy.prepare_rotamers_for_packing( pose, *rotsets.rotamer_set_for_residue( ii ) );
		for ( Size ii = 8; ii <= 11; ++ii )	{
			TS_ASSERT( rotsets.rotamer_set_for_residue(ii)->get_trie( vdw_method )() != 0 );
			//std::cout << "trie: " << ii << " ";
			//static_cast< core::scoring::trie::RotamerTrieBase const & > (*rotsets.rotamer_set_for_residue(ii)->get_trie( vdw_method ) ).print();
			//std::cout << std::endl;
		}

		Size count_comparisons( 0 ), count_wrong(0);
		for ( Size ii = 8; ii <= 11; ++ii ) {
			RotamerSet const & iiset = *rotsets.rotamer_set_for_residue( ii );
			for ( Size jj = ii+1; jj <= 11; ++jj ) {
				RotamerSet const & jjset = *rotsets.rotamer_set_for_residue( jj );
				// compute the rotamer pair energies for ii/jj
				ObjexxFCL::FArray2D< core::PackerEnergy > energy_table( jjset.num_rotamers(), iiset.num_rotamers(), core::PackerEnergy(0.0) );
				vdw_energy.evaluate_rotamer_pair_energies( iiset, jjset, pose, sfxn, sfxn.weights(), energy_table );

				/// And now verify that it worked
				ObjexxFCL::FArray2D< core::PackerEnergy > temp_table3( energy_table );
				temp_table3 = 0;
				EnergyMap emap;
				for ( Size kk = 1, kk_end = iiset.num_rotamers(); kk <= kk_end; ++kk ) {
					for ( Size ll = 1, ll_end = jjset.num_rotamers(); ll <= ll_end; ++ll ) {
						++count_comparisons;
						emap.zero();
						vdw_energy.residue_pair_energy( *iiset.rotamer( kk ), *jjset.rotamer( ll ), pose, sfxn, emap );
						temp_table3( ll, kk ) += sfxn.weights().dot( emap );
						TS_ASSERT( std::abs( energy_table( ll, kk ) - temp_table3( ll, kk )) <= 1e-3 );
						if ( std::abs( energy_table( ll, kk ) - temp_table3( ll, kk )) > 1e-3 ) {
							++count_wrong;
							std::cout << "Residues " << iiset.resid() << " & " << jjset.resid() << " rotamers: " << kk << " & " << ll;
							std::cout << " tvt/reg discrepancy: tvt= " <<  energy_table( ll, kk ) << " reg= " << temp_table3( ll, kk );
							std::cout << " delta: " << energy_table( ll, kk ) - temp_table3( ll, kk ) << std::endl;
						}
					}
				}
			}
		}
		//std::cout << "Total energy comparisons: " << count_comparisons << " total wrong " << count_wrong << std::endl;
	}

};

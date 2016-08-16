// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/FASTERInteractionGraph.cxxtest.hh
/// @brief  test suite for the FASTER interaction graph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


// Test headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// C++ headers
#include <iostream>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/AbstractRotamerTrie.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/SimpleInteractionGraph.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/pool/poolfwd.hpp>


class SimpleInteractionGraphTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	void test_instantiate_simple_ig() {

		using namespace core::chemical;
		using namespace core::graph;
		using namespace core::pack;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		using namespace core::pose;
		using namespace core::scoring;
		using core::Size;

		PoseOP trpcage = create_trpcage_ideal_poseop();
		PackerTaskOP task = TaskFactory::create_packer_task( *trpcage );
		//task->or_preserve_c_beta( true );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= 20; ++ii ) {
			if ( ii == 11 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		ScoreFunctionOP sfxn = get_score_function();
		sfxn->set_weight( fa_pair, 0.0 );
		sfxn->set_weight( hbond_sc, 0.0 );
		sfxn->set_weight( hbond_bb_sc, 0.0 );

		core::Real startscore = (*sfxn)( *trpcage ); // score the pose first;
		sfxn->setup_for_packing( *trpcage, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( *trpcage, *sfxn, task );

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( *trpcage, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( *trpcage, *sfxn );

		//trpcage->dump_pdb( "test_trpcage_1.pdb" );

		SimpleInteractionGraphOP simple_ig( new SimpleInteractionGraph );
		simple_ig->set_scorefunction( *sfxn );
		simple_ig->initialize( *trpcage );

		/*std::cout << "Initial pose one-body energy: " << sfxn->weights().dot( trpcage->energies().onebody_energies( 12 )) << std::endl;
		for ( Node::EdgeListConstIter iter = trpcage->energies().energy_graph().get_node( 12 )->edge_list_begin(),
		iter_end = trpcage->energies().energy_graph().get_node( 12 )->edge_list_end(); iter != iter_end; ++iter ) {
		std::cout << "Energy Edge to " << (*iter)->get_other_ind( 12 ) << " ";
		EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*iter);
		EnergyMap emap = eedge->fill_energy_map();
		//emap.show_weighted( std::cout, sfxn->weights() );

		std::cout << sfxn->weights().dot( emap ) << std::endl;
		}
		std::cout << "Initial one-body energy: " << static_cast< SimpleNode const * > (simple_ig->get_node( 12 ))->current_one_body_energy() << std::endl;
		for ( Node::EdgeListConstIter iter = simple_ig->get_node( 12 )->edge_list_begin(),
		iter_end = simple_ig->get_node( 12 )->edge_list_end(); iter != iter_end; ++iter ) {
		SimpleEdge const * simple_edge = static_cast< SimpleEdge const * > ( *iter );
		std::cout << "Curr energy: " << simple_edge->get_other_ind( 12 ) << " " << simple_edge->get_current_energy() << std::endl;
		}*/


		core::conformation::ResidueCOP res12rot4 = rotsets->rotamer_set_for_residue( 12 )->rotamer( 4 );

		core::Real negdelta = simple_ig->consider_substitution( 12, res12rot4 );
		trpcage->replace_residue( 12, *res12rot4, false );
		core::Real score_after_res12rot4 = (*sfxn)(*trpcage);
		//trpcage->dump_pdb( "test_trpcage_2.pdb" );

		TS_ASSERT_DELTA( negdelta, startscore - score_after_res12rot4, 1e-6 );

		/*std::cout << "Afterwards pose one-body energy: " << sfxn->weights().dot( trpcage->energies().onebody_energies( 12 )) << std::endl;
		for ( Node::EdgeListConstIter iter = trpcage->energies().energy_graph().get_node( 12 )->edge_list_begin(),
		iter_end = trpcage->energies().energy_graph().get_node( 12 )->edge_list_end(); iter != iter_end; ++iter ) {
		std::cout << "Energy Edge to " << (*iter)->get_other_ind( 12 ) << " ";
		EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*iter);
		EnergyMap emap = eedge->fill_energy_map();
		//emap.show_weighted( std::cout, sfxn->weights() );

		std::cout << sfxn->weights().dot( emap ) << std::endl;
		}
		std::cout << "Alternate one-body energy: " << static_cast< SimpleNode const * > (simple_ig->get_node( 12 ))->proposed_one_body_energy() << std::endl;
		for ( Node::EdgeListConstIter iter = simple_ig->get_node( 12 )->edge_list_begin(),
		iter_end = simple_ig->get_node( 12 )->edge_list_end(); iter != iter_end; ++iter ) {
		SimpleEdge const * simple_edge = static_cast< SimpleEdge const * > ( *iter );
		std::cout << "Alt energy: " << simple_edge->get_other_ind( 12 ) << " " << simple_edge->get_proposed_energy() << std::endl;
		}*/
	}

	void test_simple_ig_neutral_compare() {
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::conformation;
		using namespace core::pack::interaction_graph;

		core::pose::Pose  oneten(create_1ten_pdb_pose());
		ScoreFunctionOP sfxn = get_score_function();
		sfxn->score(oneten);

		SimpleInteractionGraphOP simple_ig( new SimpleInteractionGraph );
		simple_ig->set_scorefunction( *sfxn );
		simple_ig->initialize( oneten );

		for ( core::Size ii(1); ii <= oneten.n_residue(); ++ii ) {
			ResidueCOP res( ResidueOP( new Residue(oneten.residue(ii)) ) );
			TS_ASSERT_DELTA( simple_ig->consider_substitution(ii, res), 0, 0.0001);
			simple_ig->reject_change( ii );
		}
	}
};

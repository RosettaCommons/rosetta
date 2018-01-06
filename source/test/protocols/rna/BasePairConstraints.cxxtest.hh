// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/rna/BasePairConstraints.cxxtest.hh
/// @brief  test suite for BasePairConstraints -- they live in core/scoring
/// @details But I think the test should live in protocols nonetheless because
/// the way I want to write it makes heavy use of SubMotifLibrary, kinda.
/// @author Andy Watkins (amw579@stanford.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers

#include <platform/types.hh>

// Package Headers
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/BasePairConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/id/AtomID.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/leontis_westhof_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/import_pose/libraries/RNA_JumpLibrary.hh>
#include <core/import_pose/libraries/RNA_LibraryManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/kinematics/FoldTree.hh>
#include <numeric/xyzVector.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::kinematics;
using namespace core::chemical::rna;
using namespace core::pose;
using namespace core::pose::rna;
using namespace core::import_pose;
using namespace core::import_pose::libraries;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::scoring::methods;
using namespace protocols::rigid;

static basic::Tracer TR("BasePairConstraintTests");

class BasePairConstraintTests : public CxxTest::TestSuite {

public:

	void
	save_pose_as_submotif( core::pose::PoseOP const & pose, std::string const & tag ) {
		poses_[ tag ].push_back( pose );
	}

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-mute core.pose protocols.moves.RigidBodyMover core.chemical -out:level 300" );

		scorefxn_ = ScoreFunctionOP( new ScoreFunction() );
		scorefxn_->set_weight( base_pair_constraint, 1.0 );
		core::chemical::ResidueTypeSetCOP rsd_set_op = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		RNA_JumpLibraryCOP rna_jump_library( RNA_LibraryManager::get_instance()->rna_jump_library_cop() );
		for ( Size i = 0; i <= 3; i++ ) {
			for ( Size j = 0; j <= 3; j++ ) {
				for ( Size e1 = 1; e1 <= 3; e1++ ) {
					for ( Size e2 = 1; e2 <= 3; e2++ ) {
						for ( Size o = 1; o <= 2; o++ ) {
							BasePairType const base_pair_type( rna_nts[ i ], rna_nts[ j ], BaseEdge( e1 ), BaseEdge( e2 ), BaseDoubletOrientation( o ) );
							if ( !rna_jump_library->has_template( base_pair_type ) ) continue;
							//TR << "Storing base pair archetype for: " << base_pair_type.tag() << std::endl;
							RNA_PairingTemplateList const & templates( rna_jump_library->rna_pairing_template_map( base_pair_type ) );

							PoseOP pose( new Pose );
							std::string sequence;
							sequence += rna_nts[ i ];
							sequence += rna_nts[ j ];
							make_pose_from_sequence( *pose, sequence, *rsd_set_op, false /*auto termin*/  );

							FoldTree f( 2 );
							f.new_jump( 1, 2, 1 );
							fill_in_default_jump_atoms( f, *pose );
							pose->fold_tree( f );
							core::kinematics::Jump const & j( templates[ 1 ]->jump_forward() );
							pose->set_jump( 1, j );
							core::import_pose::cleanup( *pose, true /* force_cut_at_rna_chainbreak */ );

							std::string tag( base_pair_type.tag() + "_00001" );
							//pose->dump_pdb( tag + ".pdb" );
							save_pose_as_submotif( pose, tag );
						} // o
					} // e2
				} // e1
			} // j
		} // i
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_all_constraint_possibilities()
	{
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;

		for ( Size i = 0; i <= 3; i++ ) {
			for ( Size j = 0; j <= 3; j++ ) {
				for ( Size e1 = 1; e1 <= 3; e1++ ) {
					for ( Size e2 = 1; e2 <= 3; e2++ ) {
						for ( Size o = 1; o <= 2; o++ ) {
							BasePairType const base_pair_type( rna_nts[ i ], rna_nts[ j ], BaseEdge( e1 ), BaseEdge( e2 ), BaseDoubletOrientation( o ) );

							// Grab a pose with this kind of base pair.
							std::string tag( base_pair_type.tag() + "_00001" );

							if ( poses_.find( tag ) == poses_.end() ) continue;
							Pose & pose = *( poses_.at( tag )[ 1 ] );
							// Maybe it needs to be ambiguous for the test.
							/*
							* auto bpc = BasePairConstraintOP( new BasePairConstraint( 1, 2, BaseEdge( e1 ), BaseEdge( e2 ),
							* core::pose::rna::get_LW_orientation( BaseEdge( e1 ), BaseEdge( e2 ), BaseDoubletOrientation( o ) ) ) );
							*
							* bpc->init_subsidiary_constraints( pose );
							*/
							auto bpc2 = BasePairConstraintOP( new BasePairConstraint( 1, 2, BaseEdge( e1 ), BaseEdge( e2 ),
								core::pose::rna::get_LW_orientation( BaseEdge( e1 ), BaseEdge( e2 ), BaseDoubletOrientation( o ) ) ) );
							auto bpc  = BasePairConstraintOP( new BasePairConstraint( 2, 1, BaseEdge( e2 ), BaseEdge( e1 ),
								core::pose::rna::get_LW_orientation( BaseEdge( e2 ), BaseEdge( e1 ), BaseDoubletOrientation( o ) ) ) );

							auto ambig = AmbiguousConstraintOP( new AmbiguousConstraint() );
							bpc->init_subsidiary_constraints( pose );
							bpc2->init_subsidiary_constraints( pose );

							ambig->add_individual_constraint( bpc );
							ambig->add_individual_constraint( bpc2 );

							pose.add_constraint( ambig );

							Real score = ( *scorefxn_ )( pose );
							//pose.dump_pdb( tag + ".pdb" );
							//scorefxn_->show( TR, pose );
							TR << "Base pair constraint score of pose with tag " << tag << " is: " << score << std::endl;
							TS_ASSERT( score <= 1000 );
						}
					}
				}
			}
		}
	}

	void test_bad_when_stretched_all_constraint_possibilities()
	{
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		RigidBodyTransMover rbtm( core::Vector(1,1,1), 1, false );
		rbtm.step_size( 10000 );

		for ( Size i = 0; i <= 3; i++ ) {
			for ( Size j = 0; j <= 3; j++ ) {
				for ( Size e1 = 1; e1 <= 3; e1++ ) {
					for ( Size e2 = 1; e2 <= 3; e2++ ) {
						for ( Size o = 1; o <= 2; o++ ) {
							BasePairType const base_pair_type( rna_nts[ i ], rna_nts[ j ], BaseEdge( e1 ), BaseEdge( e2 ), BaseDoubletOrientation( o ) );

							// Grab a pose with this kind of base pair.
							std::string tag( base_pair_type.tag() + "_00001" );

							if ( poses_.find( tag ) == poses_.end() ) continue;
							Pose & pose = *( poses_.at( tag )[ 1 ] );
							rbtm.apply( pose );
							auto bpc = BasePairConstraintOP( new BasePairConstraint( 1, 2, BaseEdge( e1 ), BaseEdge( e2 ),
								core::pose::rna::get_LW_orientation( BaseEdge( e1 ), BaseEdge( e2 ), BaseDoubletOrientation( o ) ) ) );
							bpc->init_subsidiary_constraints( pose );
							pose.add_constraint( bpc );

							Real score = ( *scorefxn_ )( pose );
							//scorefxn_->show( TR, pose );
							//pose.dump_pdb( tag + "_stretched.pdb" );
							TR << "Base pair constraint score of pose with tag " << tag << " is: " << score << std::endl;
							TS_ASSERT( score > 1000 );
						}
					}
				}
			}
		}
	}

	core::scoring::ScoreFunctionOP scorefxn_ = nullptr;
	std::map< std::string, utility::vector1< PoseOP > > poses_;

};



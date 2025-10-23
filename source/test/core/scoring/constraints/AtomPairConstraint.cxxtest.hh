// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/AngleConstraint.cxxtest.hh
/// @brief  test suite for angle constraints
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AtomToAxisConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/FourPointsFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>


//Auto Headers
#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <sstream>
#include <string>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.AtomPairConstraint.cxxtest");

using namespace core;

class AtomPairConstraintTests : public CxxTest::TestSuite
{

public:
	AtomPairConstraintTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_atom_pair_score() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::scoring::func::FourPointsFunc fourpts;
		fourpts.xyz( 1, Vector( 0, 0, 0 ) );
		fourpts.xyz( 2, Vector( 0, 1.0, 0 ));
		fourpts.xyz( 3, Vector( 0.707, 0.707, 0 ));
		fourpts.xyz( 4, Vector( 0.707, 0.707, 1.0 )); // 90 degrees

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 1.2, 0.5 ) );

		AtomID at1( 1, 1), at2( 2, 1 );

		AtomPairConstraint atom_pair_cst( at1, at2, func );
		EnergyMap weights, emap;
		weights[ atom_pair_constraint ] = 1.0;
		atom_pair_cst.score( fourpts, weights, emap );
		//Size before_precision = std::cout.precision();
		//std::cout.precision( 16 );
		//std::cout << "Atom pair constraint func: " << emap[ atom_pair_constraint ] << std::endl;
		//std::cout.precision( before_precision );
		TS_ASSERT_DELTA( emap[ atom_pair_constraint ],   0.1599999999999999, 1e-12 );

	}

	void test_atom_pair_derivatives() {

		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		core::pose::PoseOP ubqstump = create_twores_1ubq_poseop();
		TS_ASSERT( ubqstump->size() == 2 );

		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 1.2, 0.5 ) );
		AtomID at1, at2;
		ScoreFunction sfxn;
		sfxn.set_weight( atom_pair_constraint, 1.0 );
		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		for ( Size ii = 1; ii <= ubqstump->size(); ++ii ) {
			core::chemical::ResidueType const & rsd_type( ubqstump->residue_type( ii ));
			// for each dihedral angle in the residue type
			for ( Size jj = 1; jj <= rsd_type.natoms(); ++jj ) {
				for ( Size kk = jj + 1; kk <= rsd_type.natoms(); ++kk ) {
					if ( rsd_type.path_distance( jj, kk ) != 1 ) continue;

					at1.rsd() = at2.rsd() = ii;
					at1.atomno() = jj;
					at2.atomno() = kk;

					AtomPairConstraintOP atom_pair_cst( new AtomPairConstraint( at1, at2, func ) );
					ubqstump->remove_constraints();
					ubqstump->add_constraint( atom_pair_cst );

					AtomDerivValidator adv( *ubqstump, sfxn, movemap );
					/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
					/// that the analytic norm matches the numeric norm to 1e-3.
					adv.simple_deriv_check( true, 1e-6 );
				}
			}
		}
	}


	void test_atom_pair_constraint_clone() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		AtomID at1( 1, 1), at2( 2, 1 );
		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 1.2, 0.5 ) );

		AtomPairConstraintOP atom_pair_cst( new AtomPairConstraint( at1, at2, func ));
		ConstraintOP cloned_cst = atom_pair_cst->clone();
		AtomPairConstraintOP cloned_apc = utility::pointer::dynamic_pointer_cast< AtomPairConstraint > ( cloned_cst );

		// ensure the dynamic cast succeeds
		TS_ASSERT( cloned_apc );

		// Make sure that the clone isn't the same as the original -- of course, right?
		TS_ASSERT_DIFFERS( atom_pair_cst, cloned_cst );
		TS_ASSERT_DIFFERS( atom_pair_cst, cloned_apc );

		// check mutual equality; a == b and b == a
		TS_ASSERT( *atom_pair_cst == *cloned_cst );
		TS_ASSERT( *cloned_cst == *atom_pair_cst );

		// clone() should perform a deep copy of the internal func object, verifiable by looking
		// at the func pointers and making sure they point at different objects.
		TS_ASSERT_DIFFERS( & atom_pair_cst->get_func(), & cloned_cst->get_func() );

	}

	void test_atom_pair_constraint_equality_operator() {
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		AtomID at1( 1, 1), at2( 2, 1 );
		core::scoring::func::HarmonicFuncOP func( new core::scoring::func::HarmonicFunc( 1.2, 0.5 ) );

		AtomPairConstraintOP atom_pair_cst1( new AtomPairConstraint( at1, at2, func ));

		core::scoring::func::HarmonicFuncOP func2( new core::scoring::func::HarmonicFunc( 1.2, 0.75 ) );
		AtomPairConstraintOP atom_pair_cst2( new AtomPairConstraint( at1, at2, func2 ));

		AtomID at3( 3, 1), at4( 4, 1 );
		AtomPairConstraintOP atom_pair_cst3( new AtomPairConstraint( at3, at2, func ));
		AtomPairConstraintOP atom_pair_cst4( new AtomPairConstraint( at1, at4, func ));

		// func objects differ; check mutual inequality; a != b and b != a
		TS_ASSERT( *atom_pair_cst1 != *atom_pair_cst2 );
		TS_ASSERT( *atom_pair_cst2 != *atom_pair_cst1 );

		// atom1s differ; check mutual inequality; a != b and b != a
		TS_ASSERT( *atom_pair_cst1 != *atom_pair_cst3 );
		TS_ASSERT( *atom_pair_cst3 != *atom_pair_cst1 );

		// atom2s differ; check mutual inequality; a != b and b != a
		TS_ASSERT( *atom_pair_cst1 != *atom_pair_cst4 );
		TS_ASSERT( *atom_pair_cst4 != *atom_pair_cst1 );
	}

	void test_PDBno_readin() {
		using namespace core::scoring::constraints;
		core::pose::Pose pose( create_twores_1ubq_pose() );
		pose.pdb_info()->number(1,393);
		pose.pdb_info()->chain(1,"B");
		// Yup - inverted chain order.
		pose.pdb_info()->number(2,345);
		pose.pdb_info()->chain(2,"A");

		std::stringstream infile(" CA 345A CA 393B HARMONIC 4.3 0.25 1 ");

		AtomPairConstraintOP atom_pair_cst( new AtomPairConstraint );

		atom_pair_cst->read_def( infile, pose, ConstraintIO::get_func_factory() );

		TS_ASSERT_EQUALS( atom_pair_cst->atom1().rsd(), 2);
		TS_ASSERT_EQUALS( atom_pair_cst->atom2().rsd(), 1);
	}

	void test_serialize_AtomPairConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		AtomID at3( 3, 1), at4( 4, 1 );
		AtomPairConstraintOP instance( new AtomPairConstraint( at3, at4, some_func ) );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		AtomPairConstraintOP instance2( new AtomPairConstraint() );
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

	void test_atom_to_axis(){
		using namespace core;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;

		ScoreFunction sfxn;
		TS_ASSERT_EQUALS( sfxn.get_nonzero_weighted_scoretypes().size(), 0 );
		sfxn.set_weight( atom_pair_constraint, 1.0 );
		TS_ASSERT_EQUALS( sfxn.get_nonzero_weighted_scoretypes().size(), 1 );

		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "GGG/GF", "fa_standard" );
		sfxn.score( pose );
		TS_ASSERT( pose.num_chains() == 2 );
		TS_ASSERT( pose.size() == 5 );

		core::pose::Pose pose2 = pose;

		utility::vector1< core::id::AtomID > axis1, axis2;
		axis1.emplace_back( 1, 1 );
		axis2.emplace_back( 1, 3 );
		core::id::AtomID const atom1( 1, 1 );
		core::id::AtomID const atom2( 1, 3 );

		TS_ASSERT( pose.atom_tree().has( atom1 ) );
		TS_ASSERT( pose.atom_tree().has( atom2 ) );

		//residue 4 is " CA"...?
		core::conformation::Residue const & Fres = pose.residue(5);

		//Distance to center of the ring
		core::Real const target_dist =
			Fres.xyz( "CG" ).distance( Fres.xyz( "CZ" ) ) / 2.0;

		FuncOP const func( new HarmonicFunc( target_dist, 0.0001 ) );

		for ( std::string s : { "CG", "CD1", "CD2", "CE1", "CE2", "CZ" } ) {
			core::id::AtomID const atom( Fres.atom_index(s), 5 );
			TS_ASSERT( pose.atom_tree().has( atom ) );
			AtomToAxisConstraintOP cst(new AtomToAxisConstraint( atom, axis1, axis2, func ));
			pose.add_constraint( cst );
			//                                                          v REVERSED v
			AtomToAxisConstraintOP cst2(new AtomToAxisConstraint( atom, axis2, axis1, func ));
			pose2.add_constraint( cst2 );
		}

		TS_ASSERT_EQUALS( pose.constraint_set()->get_all_constraints().size(), 6 );
		TS_ASSERT_EQUALS( pose2.constraint_set()->get_all_constraints().size(), 6 );

		TR << "Score: " << sfxn( pose ) << std::endl;

		//Axis order shouldn't matter!
		TS_ASSERT_DELTA( sfxn( pose ), sfxn( pose2 ), 0.01 );
	}

};
